#ifndef MESH_HPP
#define MESH_HPP

#include <NumCpp.hpp>

#include <vector>
#include <iterator>
#include <array>
#include <numeric>

// export
#include "splineconcepts_export.h"

#include "Spline.hpp"

/**
 * @brief A class containing a mesh of splines.
 */
class SPLINECONCEPTS_EXPORT Mesh
{
  // Construction / Destruction
public:
  Mesh() = default;
  Mesh(std::vector<std::shared_ptr<Spline>> splines) : splines_(splines) {}
  virtual ~Mesh() = default;

  // Public Methods
public:
  
  /**
   * @brief Closure check for the mesh.
   * 
   * A mesh is considered closed, when all splines are connected such that each spline has at least one perpendicular neighbor.
   * This is a simple test and does NOT check whether the mesh has holes (as it should be able to), but rather if we have dangling
   * splines entries that we cannot make a mesh grid from.
   */
  bool isClosed()
  {
    // check if all splines have at least one perpendicular neighbor
    for (size_t i = 0; i < splines_.size(); ++i)
    {
      bool hasPerp = false;
      for (size_t j = 0; j < splines_.size(); ++j)
      {
        if (i == j)
        {
          continue;
        }
        if (perpMap_[i][j])
        {
          hasPerp = true;
          break;
        }
      }
      if (!hasPerp)
      {
        return false;
      }
    }
    return true;
  }

  Mesh& clean(bool inPlace = false)
  {
    // cleaning rules: choose a direction by which the splines are split
    auto& max_intersect_spline = splines_[0];
    auto max_intersect_spline_index = 0;
    for (size_t i = 0; i < splines_.size(); ++i)
    {
      auto& spline = splines_[i];
      for (size_t j = 0; j < splines_.size(); ++j)
      {
        if (i == j)
        {
          continue;
        }
        if (perpMap_[i][j])
        {
          // check if the splines intersect
          if (spline->intersects(*splines_[j]))
          {
            // check if the intersection is larger than the current max
            if (spline->intersection(*splines_[j]).size() > max_intersect_spline->intersection(*splines_[max_intersect_spline_index]).size())
            {
              max_intersect_spline = spline;
              max_intersect_spline_index = i;
            }
          }
        }
      }
    }
    // split other splines at this spline
    for (size_t i = 0; i < splines_.size(); ++i)
    {
      if (i == max_intersect_spline_index)
      {
        continue;
      }
    }
  }

  /**
   * @brief Compute spline array hash.
   */
  std::size_t splineHash()
  {
    std::size_t hash = 0;
    for (std::shared_ptr<Spline> spline : splines_)
    {
      hash += spline->hash();
    }
    return hash;
  }

  /**
   * @brief Compute the mesh of a section.
   * A section is defined by at least two splines that are perpendicular to each other.
   */
  std::shared_ptr<Geometry> computeSection(std::vector<std::shared_ptr<Spline>> section)
  {
    // check if all splines are connected. We try to find all splines starting with the first spline only through map connection traversal.
    std::vector<std::shared_ptr<Spline>> connectedSplines;
    connectedSplines.push_back(section[0]);
    for (size_t i = 0; i < connectedSplines.size(); ++i)
    {
      for (size_t j = 0; j < section.size(); ++j)
      {
        if (std::find(connectedSplines.begin(), connectedSplines.end(), section[j]) != connectedSplines.end())
        {
          continue;
        }
        if (perpMap_[i][j])
        {
          connectedSplines.push_back(section[j]);
        }
      }
    }
    if (connectedSplines.size() != section.size())
    {
      throw std::invalid_argument("computeSection: section is not connected.");
    }
  }

private:
  std::vector<std::shared_ptr<Spline>> splines_;
  // hash of the spline object array
  std::size_t spline_hash_;
  // perpendicularity map
  std::vector<std::vector<bool>> perpMap_;
};


#endif // MESH_HPP
