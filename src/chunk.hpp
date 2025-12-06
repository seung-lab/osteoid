#ifndef __FASTOSTEOID_CHUNK_HPP__
#define __FASTOSTEOID_CHUNK_HPP__

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace fastosteoid::chunk {

template <typename T = float>
class Vec3 {
public:
  T x, y, z;
  Vec3() : x(0), y(0), z(0) {}
  Vec3(T x, T y, T z) : x(x), y(y), z(z) {}

  Vec3 operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z);
  }
  void operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
  }
  Vec3 operator+(const T other) const {
    return Vec3(x + other, y + other, z + other);
  }
  void operator+=(const T other) {
    x += other;
    y += other;
    z += other;
  }
  Vec3 operator-() const {
    return Vec3(-x,-y,-z);
  }
  Vec3 operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z);
  }
  Vec3 operator-(const T scalar) const {
    return Vec3(x - scalar, y - scalar, z - scalar);
  }
  Vec3 operator*(const T scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);
  }
  void operator*=(const T scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
  }
  Vec3 operator*(const Vec3& other) const {
    return Vec3(x * other.x, y * other.y, z * other.z);
  }
  void operator*=(const Vec3& other) {
    x *= other.x;
    y *= other.y;
    z *= other.z;
  }
  Vec3 operator/(const Vec3& other) const {
    return Vec3(x/other.x, y/other.y, z/other.z);
  }
  Vec3 operator/(const T divisor) const {
    return Vec3(x/divisor, y/divisor, z/divisor);
  }
  void operator/=(const T divisor) {
    x /= divisor;
    y /= divisor;
    z /= divisor;
  }
  bool operator==(const Vec3& other) const {
    return x == other.x && y == other.y && z == other.z;
  }
  T get(const int idx) const {
    switch (idx) {
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        throw new std::runtime_error("Index out of bounds.");
    }
  }
  T& operator[](const int idx) {
    return get(idx);
  }
  T dot(const Vec3& o) const {
    return x * o.x + y * o.y + z * o.z;
  }
  Vec3 abs() const {
    return Vec3(std::abs(x), std::abs(y), std::abs(z));
  }
  int argmax() const {
    if (x >= y) {
      return (x >= z) ? 0 : 2;
    }
    return (y >= z) ? 1 : 2;
  }
  T max() const {
    return std::max(x,std::max(y,z));
  }
  T min() const {
    return std::min(x,std::min(y,z));
  }
  float len() const {
    return sqrt(x*x + y*y + z*z);
  }
  float len2() const {
    return x*x + y*y + z*z;
  }
  Vec3 hat() const {
    const float l = len();
    Vec3 ret(x,y,z);
    if (l == 1) {
      return ret;
    }
    ret.x /= l;
    ret.y /= l;
    ret.z /= l;
    return ret;
  }
  bool close(const Vec3& o) const {
    return (*this - o).len2() < 1e-4;
  }
  Vec3 cross(const Vec3& o) const {
    return Vec3(
      y * o.z - z * o.y, 
      z * o.x - x * o.z,
      x * o.y - y * o.x
    );
  }
  bool is_null() const {
    return x == 0 && y == 0 && z == 0;
  }
  int num_zero_dims() const {
    return (x == 0) + (y == 0) + (z == 0);
  }
  int num_non_zero_dims() const {
    return (x != 0) + (y != 0) + (z != 0);
  }
  bool is_axis_aligned() const {
    return ((x != 0) + (y != 0) + (z != 0)) == 1;
  }
  void print(const std::string &name) const {
    if constexpr (std::is_same<T, float>::value) {
      printf("%s %.3f, %.3f, %.3f\n",name.c_str(), x, y, z);  
    }
    else {
      printf("%s %d, %d, %d\n",name.c_str(), x, y, z);
    }
  }
};

typedef std::pair<Vec3<float>, Vec3<float>> Line;

template <typename EDGE_T>
struct LineObject {
  std::vector<float> points;
  std::vector<EDGE_T> edges;

  void add_point(const Vec3<float>& pt) {
    points.push_back(pt.x);
    points.push_back(pt.y);
    points.push_back(pt.z);
  }

  void add_edge(
    const EDGE_T e1, 
    const EDGE_T e2
  ) {
    edges.push_back(e1);
    edges.push_back(e2);
  }

  void add_line(const Line& line) {
    EDGE_T i = last_edge();

    points.push_back(line.first.x);
    points.push_back(line.first.y);
    points.push_back(line.first.z);

    points.push_back(line.second.x);
    points.push_back(line.second.y);
    points.push_back(line.second.z);

    edges.push_back(i + 1);
    edges.push_back(i + 2);
  }

  // Note: exploits underflow and +1 in add_line
  EDGE_T last_edge() const {
    return (points.size() > 0) 
      ? ((points.size() - 1) / 3)
      : -1;
  }
};


Vec3<int64_t> zone2grid(int64_t zone, const Vec3<int64_t>& gs) {
    int64_t z = zone / gs.x / gs.y;
    int64_t y = (zone - gs.x * gs.y * z) / gs.x;
    int64_t x = zone - gs.x * (y + gs.y * z);
    return Vec3<int64_t>(x,y,z);
}

Vec3<float> intersect(int axis, float plane_offset, const Vec3<float> &p, const Vec3<float> &q) {
  float t = (plane_offset - p.get(axis)) / (q.get(axis) - p.get(axis));
  return p + (q - p) * t;
}

std::vector<Line> divide_line(
  const int axis,
  const float plane_value,
  const Vec3<float>& v1,
  const Vec3<float>& v2
) {
    uint8_t above[2] = {};
    uint8_t below[2] = {};
    uint8_t equal[2] = {};
    const Vec3<float> verts[2] = { v1, v2 };

    constexpr float epsilon = 1e-5;
    int aboveCount = 0, belowCount = 0, equalCount = 0;

    uint8_t i = 0;
    for (const auto& vertex : verts) {
      const float dist = vertex.get(axis) - plane_value;

      if (std::abs(dist) < epsilon) {
        equal[equalCount++] = i;
      }
      else if (dist > 0) {
        above[aboveCount++] = i;
      } 
      else {
        below[belowCount++] = i;
      }
      i++;
    }

    std::vector<Line> result;

    if (!(aboveCount == 1 && belowCount == 1)) {
      result.emplace_back(v1, v2);
    }
    else {
      const Vec3<float>& a = verts[above[0]];
      const Vec3<float>& b = verts[below[0]];

      Vec3<float> i1 = intersect(axis, plane_value, a, b);;
      result.emplace_back(a, i1);
      result.emplace_back(i1, b);
    }

    return result;
}

std::vector<Line> divide_line(
  const int axis,
  const float plane_value,
  const Line& line
) {
  return divide_line(axis, plane_value, line.first, line.second);
}

// more elegant algorithmically, but not the fastest or simpliest
// division of the triangle into subtriangles
template <typename EDGE_T>
void resect_line_iterative(
  const float* vertices,
  const Vec3<float>& minpt,
  const std::vector<int64_t>& zones,
  std::vector<LineObject<EDGE_T>>& line_grid,
  const Vec3<float>& cs,
  const Vec3<int64_t>& gs,
  const uint64_t e1,
  const uint64_t e2
) {
  const Vec3 v1(vertices[3*e1+0], vertices[3*e1+1], vertices[3*e1+2]);
  const Vec3 v2(vertices[3*e2+0], vertices[3*e2+1], vertices[3*e2+2]);

  auto z1 = zones[e1];
  auto z2 = zones[e2];

  Vec3<int64_t> g1 = zone2grid(z1, gs);
  Vec3<int64_t> g2 = zone2grid(z2, gs);

  uint64_t gxs = std::min(g1.x, g2.x);
  uint64_t gxe = std::max(g1.x, g2.x);

  uint64_t gys = std::min(g1.y, g2.y);
  uint64_t gye = std::max(g1.y, g2.y);

  uint64_t gzs = std::min(g1.z, g2.z);
  uint64_t gze = std::max(g1.z, g2.z);

  std::vector<Line> cur_lines;
  std::vector<Line> next_lines;
  
  cur_lines.emplace_back(v1, v2);

  for (uint64_t x = gxs; x <= gxe; x++) {
    float xplane = minpt.x + x * cs.x;
    for (const auto& line : cur_lines) {
      auto new_lines = divide_line(0, xplane, line);
      next_lines.insert(next_lines.end(), new_lines.begin(), new_lines.end());
    }
    std::swap(cur_lines, next_lines);
    next_lines.clear();
  }

  for (uint64_t y = gys; y <= gye; y++) {
    float yplane = minpt.y + y * cs.y;
    for (const auto& line : cur_lines) {
      auto new_lines = divide_line(1, yplane, line);
      next_lines.insert(next_lines.end(), new_lines.begin(), new_lines.end());
    }
    std::swap(cur_lines, next_lines);
    next_lines.clear();
  }

  for (uint64_t z = gzs; z <= gze; z++) {
    float zplane = minpt.z + z * cs.z;
    for (const auto& line : cur_lines) {
      auto new_lines = divide_line(2, zplane, line);
      next_lines.insert(next_lines.end(), new_lines.begin(), new_lines.end());
    }
    std::swap(cur_lines, next_lines);
    next_lines.clear();
  }

  const float icx = 1 / cs.x;
  const float icy = 1 / cs.y;
  const float icz = 1 / cs.z;

  auto zonefn = [&](const Vec3<float>& pt) {
    uint64_t ix = static_cast<uint64_t>((pt.x - minpt.x) * icx);
    uint64_t iy = static_cast<uint64_t>((pt.y - minpt.y) * icy);
    uint64_t iz = static_cast<uint64_t>((pt.z - minpt.z) * icz);

    ix = std::min(std::max(ix, static_cast<uint64_t>(0)), static_cast<uint64_t>(gs.x - 1));
    iy = std::min(std::max(iy, static_cast<uint64_t>(0)), static_cast<uint64_t>(gs.y - 1));
    iz = std::min(std::max(iz, static_cast<uint64_t>(0)), static_cast<uint64_t>(gs.z - 1));

    return ix + gs.x * (iy + gs.y * iz);
  };

  for (const auto& line : cur_lines) {
    // v1 guaranteed to not be a border point (unless the triangle is degenerate)
    uint64_t z1 = zonefn(line.first);
    uint64_t z2 = zonefn(line.second);

    uint64_t zone = std::min(z1,z2);

    line_grid[zone].add_line(line);
  }
}

// cx = chunk size x, etc
template <typename EDGE_T>
std::vector<LineObject<EDGE_T>> chunk_line(
  const float* vertices, 
  const uint64_t num_vertices,
  const EDGE_T* edges,
  const uint64_t num_edges,
  const float cx, const float cy, const float cz,
  const std::optional<float> origin_x = std::nullopt, 
  const std::optional<float> origin_y = std::nullopt, 
  const std::optional<float> origin_z = std::nullopt
) {

  if (cx <= 0 || cy <= 0 || cz <= 0) {
    throw std::runtime_error("Chunk size must have a positive non-zero volume.");
  }

  const Vec3 cs(cx,cy,cz);

  float min_x = INFINITY;
  float min_y = INFINITY;
  float min_z = INFINITY;
  float max_x = -INFINITY;
  float max_y = -INFINITY;
  float max_z = -INFINITY;

  for (uint64_t i = 0; i < num_vertices * 3; i += 3) {
    min_x = std::min(min_x, vertices[i]);
    max_x = std::max(max_x, vertices[i]);

    min_y = std::min(min_y, vertices[i+1]);
    max_y = std::max(max_y, vertices[i+1]);

    min_z = std::min(min_z, vertices[i+2]);
    max_z = std::max(max_z, vertices[i+2]);
  }

  if (origin_x.has_value()) {
    min_x = *origin_x;
  }
  if (origin_y.has_value()) {
    min_y = *origin_y;
  }
  if (origin_z.has_value()) {
    min_z = *origin_z;
  }

  const Vec3<float> minpt(min_x, min_y, min_z);

  const int64_t gx = std::max(static_cast<int64_t>(std::ceil((max_x - min_x) / cx)), static_cast<int64_t>(1));
  const int64_t gy = std::max(static_cast<int64_t>(std::ceil((max_y - min_y) / cy)), static_cast<int64_t>(1));
  const int64_t gz = std::max(static_cast<int64_t>(std::ceil((max_z - min_z) / cz)), static_cast<int64_t>(1));

  const Vec3<int64_t> gs(gx,gy,gz);

  std::vector<int64_t> zones(num_vertices);

  const float icx = 1 / cx;
  const float icy = 1 / cy;
  const float icz = 1 / cz;

  for (uint64_t i = 0, j = 0; j < num_vertices; i += 3, j++) {
    int64_t ix = static_cast<int64_t>((vertices[i] - min_x) * icx) ;
    int64_t iy = static_cast<int64_t>((vertices[i+1] - min_y) * icy);
    int64_t iz = static_cast<int64_t>((vertices[i+2] - min_z) * icz);

    ix = std::min(std::max(ix, static_cast<int64_t>(0)), static_cast<int64_t>(gx - 1));
    iy = std::min(std::max(iy, static_cast<int64_t>(0)), static_cast<int64_t>(gy - 1));
    iz = std::min(std::max(iz, static_cast<int64_t>(0)), static_cast<int64_t>(gz - 1));

    zones[j] = ix + gx * (iy + gy * iz);
  }

  std::vector<LineObject<EDGE_T>> line_grid(gx * gy * gz);
  
  for (uint64_t i = 0; i < num_edges * 2; i += 2) {
    auto e1 = edges[i+0];
    auto e2 = edges[i+1];

    resect_line_iterative(
      vertices, minpt, 
      zones,
      line_grid, cs, gs,
      e1, e2
    );
  }

  return line_grid;
}


};

#endif
