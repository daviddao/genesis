#ifndef GENESIS_PLACEMENT_SIMULATOR_H_
#define GENESIS_PLACEMENT_SIMULATOR_H_

/**
 * @brief
 *
 * @file
 * @ingroup placement
 */

#include <stddef.h>

namespace genesis {

// =============================================================================
//     Forward Declarations
// =============================================================================

class PlacementMap;

// =============================================================================
//     Placement Simulator
// =============================================================================


/**
 * @brief Simulates Placements on a Tree.
 */
class PlacementSimulator
{
public:

    static void Generate          (PlacementMap placements, size_t n);
    static void GenerateInSubtree (PlacementMap placements, size_t n);

};

} // namespace genesis

#endif // include guard
