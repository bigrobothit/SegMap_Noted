#ifndef SEGMATCH_IMPL_DYNAMIC_VOXEL_GRID_HPP_
#define SEGMATCH_IMPL_DYNAMIC_VOXEL_GRID_HPP_

#include "segmatch/dynamic_voxel_grid.hpp"

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>

#include <glog/logging.h>
#include <laser_slam/benchmarker.hpp>
#include <pcl/common/centroid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PointIndices.h>

#include "segmatch/common.hpp"

namespace segmatch
{

// Force the compiler to reuse instantiations provided in dynamic_voxel_grid.cpp
extern template class DynamicVoxelGrid<PclPoint, MapPoint>;

//=================================================================================================
//    DynamicVoxelGrid public methods implementation
//=================================================================================================

/// \brief:  Insert Point Cloud
/// \param new_cloud: 待插入点云
/// \return: 点云中每个点在DVG中的voxel indice vector
template <_DVG_TEMPLATE_DECL_>
std::vector<int> DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::insert(const InputCloud &new_cloud)
{
    // 新点云创建对应点个数的下标数组
    std::vector<int> created_voxel_indices;
    if (new_cloud.empty())
        return created_voxel_indices;

    // 内存预留
    created_voxel_indices.reserve(new_cloud.size());

    // 点云排序，得到排序后的索引点云，此时给点云中的每个点都找到了应该属于的voxel index
    IndexedPoints_ new_points = indexAndSortPoints_(new_cloud);

    // Create containers and reserve space to prevent reallocation
    std::vector<Voxel_> new_voxels; // 加入点云之后的新的体素，按照最大需求进行内存预留，假设每个点都是一个单独的voxel
    new_voxels.reserve(voxels_.size() + new_cloud.size());
    std::unique_ptr<VoxelCloud> new_active_centroids(new VoxelCloud()); // 激活状态点中心，内存预留
    new_active_centroids->reserve(active_centroids_->size() + new_cloud.size());
    std::unique_ptr<VoxelCloud> new_inactive_centroids(new VoxelCloud()); //非激活状态点中心，内存预留
    new_inactive_centroids->reserve(inactive_centroids_->size() + new_cloud.size());

    // 设置索引点云new_points和voxels的迭代器
    auto p_it = new_points.begin();
    auto v_it = voxels_.begin();
    const auto p_end = new_points.end();
    const auto v_end = voxels_.end();

    // Merge points updating the affected voxels.
    // 终止条件是点云和voxel都遍历完毕
    while (!(p_it == p_end && v_it == v_end))
    {
        // 遍历过程中voxel data和voxel index记录
        VoxelData_ voxel_data = {nullptr, p_it, p_it}; //存放voxel_data起始点在索引点云中的位置points_begin
        IndexT voxel_index;

        // Use the next voxel if it has the upcoming index.
        // 点云遍历完毕，接着遍历voxel中的点
        // 或者当前voxel的index小于点云点需要的voxel index，也就是还没有找到当前点云应该在的那个voxel，接着遍历voxel中的点
        if ((p_it == p_end) || (v_it != v_end && v_it->index <= p_it->voxel_index))
        {
            voxel_index = v_it->index; // 当前voxel index记录
            voxel_data.old_voxel = &(*v_it);
            ++v_it; // 指向下一个voxel
        }
        else // 点云没有遍历完毕，并且找到了当前点应该属于的voxel index
        {
            voxel_index = p_it->voxel_index;
        }

        // Gather all the points that belong to the current voxel
        //
        while (p_it != p_end && p_it->voxel_index == voxel_index)
        {
            ++p_it;
        }
        voxel_data.points_end = p_it; // 存储voxel_data结束点在索引点云中的位置points_end

        // Create the voxel
        if (createVoxel_(voxel_index, voxel_data, new_voxels,
                         *new_active_centroids, *new_inactive_centroids))
        {
            created_voxel_indices.push_back(new_active_centroids->size() - 1); // 由inactive变成active的voxel的在active_centroids_下的index
        }
    }

    // Done! Save the new voxels and return the indices of the triggered
    // voxels.
    voxels_ = std::move(new_voxels);
    active_centroids_ = std::move(new_active_centroids);
    inactive_centroids_ = std::move(new_inactive_centroids);
    return created_voxel_indices;
}

/// \brief: 根据点云坐标点计算每一个点的index
template <_DVG_TEMPLATE_DECL_>
template <typename PointXYZ_>
inline IndexT DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::getIndexOf(const PointXYZ_ &point) const
{
    // 确认每个点都具有xyz三个坐标值
    static_assert(pcl::traits::has_xyz<PointXYZ_>::value,
                  "PointXYZ_ must be a structure containing XYZ coordinates");
    // TODO: One could pack indexing transformation, offsetting and scaling in a single
    // transformation. min_corner and max_corner would need to be transformed as well in order
    // to allow checks. Since it would decrease readability significantly, this should be done only
    // if optimization is really necessary.

    // Transform the point back to the grid frame for hashing.
    // 将点云从sensor frame变换到map frame
    Eigen::Vector3f transformed_coords = indexing_transformation_.transform(point.getVector3fMap());

    // Ensure that the transformed point lies inside the grid.
    // 确认每个点都在珊格范围内
    CHECK(min_corner_(0) <= transformed_coords.x() && transformed_coords.x() < max_corner_(0));
    CHECK(min_corner_(1) <= transformed_coords.y() && transformed_coords.y() < max_corner_(1));
    CHECK(min_corner_(2) <= transformed_coords.z() && transformed_coords.z() < max_corner_(2));

    // Compute voxel index of the point.
    // 转换到珊格坐标系下
    Eigen::Vector3f grid_coords = (transformed_coords + indexing_offset_) * world_to_grid_;
    // 根据论文公式（1）计算珊格坐标系下标索引
    return static_cast<IndexT>(grid_coords[0]) + (static_cast<IndexT>(grid_coords[1]) << bits_x) + (static_cast<IndexT>(grid_coords[2]) << (bits_x + bits_y));
}

template <_DVG_TEMPLATE_DECL_>
void DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::transform(
    const kindr::minimal::QuatTransformationTemplate<float> &transformation)
{
    BENCHMARK_BLOCK("SM.TransformLocalMap.TransformDVG");

    // Update transforms
    pose_transformation_ = transformation * pose_transformation_;
    indexing_transformation_ = pose_transformation_.inverse();

    // Transform point clouds in-place
    for (auto centroids : {std::ref(active_centroids_), std::ref(inactive_centroids_)})
    {
        for (auto &point : *centroids.get())
        {
            point.getVector3fMap() = transformation.transform(point.getVector3fMap());
        }
    }
}

template <_DVG_TEMPLATE_DECL_>
void DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::clear()
{
    // Reset transformations.
    pose_transformation_.setIdentity();
    indexing_transformation_.setIdentity();

    // Clear points and voxels.
    active_centroids_->clear();
    inactive_centroids_->clear();
    voxels_.clear();
}

template <_DVG_TEMPLATE_DECL_>
void DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::dumpVoxels() const
{
    for (const Voxel_ &v : voxels_)
    {
        LOG(INFO) << "Voxel " << uint32_t(v.index) << ": " << v.num_points << " " << *(v.centroid);
    }
}

//=================================================================================================
//    DynamicVoxelGrid private methods implementation
//=================================================================================================

/// \brief: 给点云计算索引，并将点云按照voxel index升序进行排序后返回IndexedPoints_形式索引点云
template <_DVG_TEMPLATE_DECL_>
inline typename DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::IndexedPoints_
DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::indexAndSortPoints_(const InputCloud &points) const
{
    IndexedPoints_ indexed_points; // vector<IndexedPoint_>，IndexedPoint_包含point和voxel_index
    indexed_points.reserve(points.size());

    // 将点云中每个点计算索引后插入到IndexedPoints_中，使用vector<IndexPoint_>存储
    for (const auto &point : points)
    {
        indexed_points.emplace_back(point, getIndexOf(point));
    }

    // 定义按照索引比较的lambda函数，并按照voxel index升序进行排序
    auto predicate = [](const IndexedPoint_ &a, const IndexedPoint_ &b) {
        return a.voxel_index < b.voxel_index;
    };
    std::sort(indexed_points.begin(), indexed_points.end(), predicate);

    return indexed_points;
}

template <_DVG_TEMPLATE_DECL_>
inline bool DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::createVoxel_(
    const IndexT index, const VoxelData_ &data,
    std::vector<Voxel_> &new_voxels, VoxelCloud &new_active_centroids,
    VoxelCloud &new_inactive_centroids)
{
    VoxelPointT centroid;
    auto centroid_map = centroid.getVector3fMap(); //使用的Eigen::Map操作，相当于引用
    uint32_t old_points_count = 0u;
    uint32_t new_points_count = std::distance(data.points_begin, data.points_end);

    // Add contribution from the existing voxel.
    if (data.old_voxel != nullptr)
    {
        centroid = *(data.old_voxel->centroid);
        old_points_count = data.old_voxel->num_points;
        if (new_points_count != 0u)
        {
            centroid_map *= static_cast<float>(old_points_count);
        }
    }
    uint32_t total_points_count = old_points_count + new_points_count;

    // Add contribution from the new points.
    if (new_points_count != 0u)
    {
        for (auto it = data.points_begin; it != data.points_end; ++it)
        {
            centroid_map += it->point.getVector3fMap();
        }
        centroid_map /= static_cast<float>(total_points_count);
    }

    // Save centroid to the correct point cloud.
    VoxelPointT *centroid_pointer;
    bool is_new_voxel = false;
    if (total_points_count >= min_points_per_voxel_) // 当前是active
    {
        new_active_centroids.push_back(centroid);
        centroid_pointer = &new_active_centroids.back();
        is_new_voxel = (old_points_count < min_points_per_voxel_); // 是否是一个新增加的voxel
    }
    else
    {
        new_inactive_centroids.push_back(centroid);
        centroid_pointer = &new_inactive_centroids.back();
    }

    new_voxels.emplace_back(centroid_pointer, index, total_points_count);
    return is_new_voxel;
}

template <_DVG_TEMPLATE_DECL_>
inline std::vector<bool> DynamicVoxelGrid<_DVG_TEMPLATE_SPEC_>::removeCentroids_(
    VoxelCloud &target_cloud, std::vector<VoxelPointT *> to_remove)
{
    std::vector<bool> is_removed(target_cloud.size(), false);
    if (to_remove.empty())
        return is_removed;

    size_t next_removal_index = 0u;
    size_t centroid_index = 0u;

    // Push one more element so that we don't read past the end of the vector.
    to_remove.push_back((VoxelPointT *)0);

    // Remove the required centroids and keep track of their indices.
    auto new_end = std::remove_if(target_cloud.begin(), target_cloud.end(),
                                  [&](VoxelPointT &p) {
                                      const bool remove_p = (&p == to_remove[next_removal_index]);
                                      if (remove_p)
                                      {
                                          is_removed[centroid_index] = true;
                                          ++next_removal_index;
                                      }
                                      ++centroid_index;
                                      return remove_p;
                                  });
    target_cloud.erase(new_end, target_cloud.end());

    return is_removed;
}

} // namespace segmatch

#endif // SEGMATCH_IMPL_DYNAMIC_VOXEL_GRID_HPP_
