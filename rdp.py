"""
This is a Ramer-Douglas-Peucker implementation that only expands the shape.
In other words, it guarantees that the original polygon is strictly contained within the new (simplified) polygon.
It does this by a slight modification of RDP: whenever we remove a point, if that point is outside the polygon,
shift the line parallel to itself to intersect that point.

There are no guarantees on whether the points of the simplified polygon are a subset of the original set of points.
Indeed, in many cases they won't be.

Usage:
  rdp(points, epsilon)
    - points should be a 2d numpy array with dimensions (n, 2). Moreover, points[-1] should equal points[0] -- i.e. points should form a closed polygon.
    - epsilon is the epsilon parameter in RDP
    - The return value is a 2d numpy array that approximates the original polygon while containing fewer segments. It may not be closed.

  rdp_closed(points, epsilon)
    - Just like rdp but guarantees that the return value will also be a closed polygon.
      This is most likely what you want most of the time.
"""


import numpy as np


def line_dists(points, start, end):
  """
  Returns the signed distances of each point to start and end
  """
  if np.all(start == end):
    return np.linalg.norm(points - start, axis=1)

  vec = end - start
  cross = np.cross(vec, start - points)
  return np.divide(cross, np.linalg.norm(vec))


def glue(seg1, seg2):
  """
  Glues two segments together
  In other words, given segments A and B which have endpoints
  that are close together, computes a "glue point" that both segments
  can be extended to in order to intersect.
  Assumes seg1 and seg2 are arrays containing two points each,
  and that if seg1 = [a, b], it can be extended in the direction of b,
  and if seg2 = [c, d], it can be extended in the direction of c.
  """

  x1 = seg1[0]
  dir1 = seg1[1] - x1
  x2 = seg2[0]
  dir2 = seg2[1] - x2

  # The set of points that are on segment 1 treated as a ray is x1 + dir1 * t for all t >= 0
  # Similarly for segment 2, but for s <= 1
  # We set up the system:
  #     x1 + dir1 * t = x2 + dir2 * s
  # --> dir1 * t + dir2 * (-s) = x2 - x1
  # --> [dir1, dir2] * [t, -s] = x2 - x1
  mat = np.matrix([dir1, dir2]).T
  # Solve for the vector [t, -s]
  try:
    t_s = np.array((np.linalg.inv(mat) @ (x2 - x1))).flatten()
    # Recall that we can't make a segment go in a backwards direction
    # So we must have t >= 0 and s <= 1. However, since we solved for [t, -s]
    # we want that t_s[0] >= 0 and t_s[1] >= -1. If this fails, set t_s to None
    if (t_s[0] < 0) or (t_s[1] < -1):
      t_s = None
  except np.linalg.LinAlgError:
    # Singular matrix i.e. parallel
    t_s = None
  if t_s is None:
    # Just connect them with a line
    return np.array([seg1[1], seg2[0]])
  else:
    # x1 + dir1 * t_s[0] and x2 - dir2 * t_s[1] Should be the same
    return np.array([x1 + dir1 * t_s[0]])


def __rdp(points, epsilon):
  """Computes expansion only rdp assuming a clockwise orientation"""
  start, end = points[0], points[-1]
  dists = line_dists(points, start, end)

  # First compute the largest point away from the line just like the ordinary RDP
  index = np.argmax(np.abs(dists))
  dmax = abs(dists[index])

  if dmax > epsilon:
    result1 = __rdp(points[:index + 1], epsilon)
    result2 = __rdp(points[index:], epsilon)
    result = np.vstack((result1[:-1], glue(result1[-2:], result2[:2]), result2[1:]))
  else:
    # All points are within epsilon of the line
    # We take the point furthest *outside* the boundary (i.e. with negative distance)
    # and shift the line segment parallel to itself to intersect that point
    new_start, new_end = np.copy(start), np.copy(end)
    vec_x, vec_y = (end - start)
    norm = np.sqrt(vec_x * vec_x + vec_y * vec_y)
    if norm != 0:
      vec_rot90 = np.array([-vec_y, vec_x]) / norm
      # TODO -- verify that this works: can optimize so that if dists[index] < 0, no need to search again, we have index_min = index and dmin = -dmax
      index_min = np.argmin(dists)
      dmin = -dists[index_min]
      if dmin > 0:
        vec_rot90 *= dmin
        new_start += vec_rot90
        new_end += vec_rot90
    result = np.array([new_start, new_end])

  return result


def rdp(points, epsilon=0):
  """Computes expansion only rdp with no assumptions on orientation"""
  if orientation(points):
    return __rdp(points, epsilon)
  else:
    return __rdp(points[::-1], epsilon)[::-1]


def rdp_closed(points, epsilon=0):
  """Ensure that expansion only rdp returns a closed loop at the end"""
  new_points = rdp(points, epsilon)
  glue_pts = glue(new_points[-2:], new_points[:2])
  if len(glue_pts) == 1:
    return np.vstack((glue_pts, new_points[1:-1], glue_pts))
  else:
    return np.vstack((new_points, [new_points[0]]))


def orientation(points):
  """Returns true if the points area in a clockwise order"""
  x, y = points.T
  # Compute signed area and compare with 0
  return np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)) >= 0

def test():
  import matplotlib.pyplot as plt
  def make_points():
    """Returns a crescent shaped thing"""
    dx = 0.03
    xs_1 = np.arange(-1, 0, dx)
    xs_2 = np.arange(0, -0.5, -dx)
    points_1 = np.array([xs_1, np.sqrt(1-xs_1*xs_1)]).T
    points_2 = np.array([xs_2, np.sqrt(1.25**2 - (xs_2-0.75)**2)]).T
    points_3 = (points_2 * np.array([1, -1]))[::-1]
    points_4 = (points_1 * np.array([1, -1]))[::-1]
    points = np.vstack((points_1, points_2, points_3, points_4))
    return points

  def pprint_numpy(points):
    for x, y in points:
      print(f'({x:.2f}, {y:.2f})', end=', ')
    print("\n")

  def plot_points(pts):
    plt.plot(pts.T[0], pts.T[1])


  pts = make_points()[::-1]
  plot_points(pts)
  pts_rdp = rdp_closed(pts, 0.2)
  plot_points(pts_rdp)
  print("Original:")
  pprint_numpy(pts)
  print("\nReduced:")
  pprint_numpy(pts_rdp)

  plt.show()

if __name__ == "__main__":
  test()

