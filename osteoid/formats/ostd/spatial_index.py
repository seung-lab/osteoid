from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

@dataclass
class OstdSpatialIndex:
  minpt:npt.NDArray[np.float32]
  maxpt:npt.NDArray[np.float32]
  chunk_size:npt.NDArray[np.float32]
  num_paths:npt.NDArray[np.uint32]
  paths:npt.NDArray[np.uint32]

  def grid_size(self) -> npt.NDArray[np.uint64]:
    grid_size = (self.maxpt - self.minpt) / self.chunk_size
    return np.ceil(grid_size).astype(np.uint64)

  def point_to_grid(self, x:float, y:float, z:float) -> int:
    grid_pt = ((pt - self.minpt) / self.chunk_size).round()
    gshape = self.grid_size()
    return grid_pt[0] + gshape[0] * (grid_pt[1] + gshape[1] * grid_pt[2])

  def ball_query(self, pt:np.ndarray, r:float):

    offsets = np.cumsum(self.num_paths)

    grid_min = self.point_to_grid(pt[0] - r, pt[1] - r, pt[2] - r)
    grid_max = self.point_to_grid(pt[0] + r, pt[1] + r, pt[2] + r)

    for gz in range(grid_min[2], grid_max[2] + 1):
      for gy in range(grid_min[1], grid_max[1] + 1):
        for gx in range(grid_min[0], grid_max[0] + 1):
          pass






