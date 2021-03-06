/**
 * @brief Implementation of TreeNode functions.
 *
 * For reasons of readability, in this implementation file, the template data types
 * NodeDataType and EdgeDataType are abbreviated NDT and EDT, respectively.
 *
 * @file
 * @ingroup tree
 */

// =============================================================================
//     Inline definitions with dependecies
// =============================================================================

// The following are inline definitions that would create circular dependecies when included in the
// class definition. Thus, they need to be here, after the definition, so that their dependend
// source files can be included without circles.
// See http://www.cplusplus.com/forum/articles/10627/ for more information on this.

// Usually, this part would be at the end of the header file, but as this is a class template,
// where the implementation (this file) is included from within the header anyway, we can as well
// outsource the inline functions to here.

#include "tree/tree_link.hpp"

namespace genesis {

/**
 * @brief True iff the node is a leaf/tip.
 */
template <class NDT, class EDT>
inline bool TreeNode<NDT, EDT>::IsLeaf() const
{
    return link_->IsLeaf();
}

/**
 * @brief True iff the node is an inner node.
 */
template <class NDT, class EDT>
inline bool TreeNode<NDT, EDT>::IsInner() const
{
    return link_->IsInner();
}

} // namespace genesis

// =============================================================================
//     Member Functions
// =============================================================================

#include "tree/tree_link.hpp"

namespace genesis {

/**
 * @brief Rank of the node, i.e. how many immediate children it has.
 */
template <class NDT, class EDT>
int TreeNode<NDT, EDT>::Rank() const
{
    int rank = -1;
    TreeLink<NDT, EDT>* link = link_;

    do {
        ++rank;
        link = link->Next();
    } while (link != link_);

    return rank;
}

/**
 * @brief Returns a one-line dump summary of the data of this node.
 */
template <class NDT, class EDT>
std::string TreeNode<NDT, EDT>::Dump() const
{
    return "Rank: " + std::to_string(Rank()) + " \t " + "name: " + this->name;
}

} // namespace genesis
