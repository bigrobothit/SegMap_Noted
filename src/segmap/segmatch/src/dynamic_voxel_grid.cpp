#include "segmatch/impl/dynamic_voxel_grid.hpp"

#include "segmatch/common.hpp"

namespace segmatch {
// Instantiate DynamicVoxelGrid for the template parameters used in the application.
template class DynamicVoxelGrid<PclPoint, MapPoint>;
pcl::PointXYZ p;
p.getVector3fMap();
// Add any other required instantiation here or in a separate file and declare them in
// segmatch/impl/dynamic_voxel_grid.hpp.
} // namespace segmatch
