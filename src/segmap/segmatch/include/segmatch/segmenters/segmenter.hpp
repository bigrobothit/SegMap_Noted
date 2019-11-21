#ifndef SEGMATCH_SEGMENTER_HPP_
#define SEGMATCH_SEGMENTER_HPP_

#include <vector>

#include "segmatch/common.hpp"
#include "segmatch/points_neighbors_providers/points_neighbors_provider.hpp"

namespace segmatch
{

// Forward declaration to speed up compilation time.
class SegmentedCloud;

/// \brief Interface for point cloud segmenters in SegMatch..
template <typename ClusteredPointT>
class Segmenter
{
  public:
    typedef pcl::PointCloud<ClusteredPointT> ClusteredCloud;

    /// \brief Finalizes an instance of the Segmenter class.
    virtual ~Segmenter() = default;

    /// \brief Cluster the given point cloud, writing the found segments in the segmented cloud.
    ///        根据所给点云，进行segmentation，并写入segmented cloud
    ///        如果cluster IDs改变了， cluster_ids_to_segment_ids 的映射也会更新
    /// If cluster IDs change, the \c cluster_ids_to_segment_ids mapping is updated accordingly.
    /// -----------------------------------------------------------------------------
    /// \param normals The normal vectors of the point cloud. This can be an empty cloud if the
    /// the segmenter doesn't require normals.
    /// \param is_point_modified Indicates for each point if it has been modified such that its
    /// cluster assignment may change.
    /// \param cloud The point cloud that must be segmented.
    /// \param points_neighbors_provider Object providing nearest neighbors information.
    /// \param segmented_cloud Cloud to which the valid segments will be added.
    /// \param cluster_ids_to_segment_ids Mapping between cluster IDs and segment IDs. Cluster
    /// \c i generates segment \c cluster_ids_to_segments_ids[i]. If
    /// \c cluster_ids_to_segments_ids[i] is equal to zero, then the cluster does not contain enough
    /// points to be considered a segment.
    /// \param renamed_segments Vectors containing segments that got a new ID, e.g. after merging
    /// two or more segments. The ordering of the vector represents the sequence of renaming
    /// operations. The first ID in each pair is the renamed segments, the second ID is the new
    /// segment ID.
    virtual void segment(
        // 点云的法线
        const PointNormals &normals,
        // 指示每个点是否已被修改以使其cluster assignment可能更改
        const std::vector<bool> &is_point_modified,
        // 待分割的点云
        ClusteredCloud &cloud,
        // 提供点云最近邻信息
        PointsNeighborsProvider<ClusteredPointT> &points_neighbors_provider,
        // 添加有效segment的点云
        SegmentedCloud &segmented_cloud,
        // cluster_id=i 和 segment_id=cluster_ids_to_segment_ids[i]的一个映射
        // 如果segment_id=0，表示该cluster包含的点不足以构成一个segment
        std::vector<Id> &cluster_ids_to_segment_ids,
        // 得到一个新segment ID的时候，包含新的segments
        // 在merge几个segments之后，vector存储了rename的顺序
        // 第一个id是renamed的segment ID，第二个ID是新的segment ID
        std::vector<std::pair<Id, Id>> &renamed_segments) = 0;
}; // class Segmenter

} // namespace segmatch

#endif // SEGMATCH_SEGMENTER_HPP_
