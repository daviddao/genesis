#ifndef GENESIS_TREE_TREENODE_H_
#define GENESIS_TREE_TREENODE_H_

/**
 * @brief This class represents a node inside of a tree.
 *
 * For more information, see TreeNode class.
 *
 * @file
 * @ingroup tree
 */

#include <string>

#include "tree/newick_processor.hpp"

namespace genesis {

// =============================================================================
//     Forward declarations
// =============================================================================

template <class NodeDataType, class EdgeDataType>
class  Tree;

template <class NodeDataType, class EdgeDataType>
class  TreeLink;

template <typename LinkPointerType, typename NodePointerType, typename EdgePointerType>
class TreeNodeIteratorLinks;

// =============================================================================
//     DefaultNodeData
// =============================================================================

class DefaultNodeData
{
public:
    /**
     * Name of the node. In case it is a leaf, this is usually the name of
     * the taxon represented by the node.
     */
    std::string name;

    /**
     * @brief Fills the node with data from a NewickBrokerElement.
     */
    inline void FromNewickBrokerElement (NewickBrokerElement* node)
    {
        name = node->name;
    }

    inline void ToNewickBrokerElement (NewickBrokerElement* node) const
    {
        node->name = name;
    }

    inline std::string Dump() const
    {
        return "Name: '" + name + "'";
    }
};

// =============================================================================
//     TreeNode
// =============================================================================

template <class NodeDataType, class EdgeDataType>
class TreeNode : public NodeDataType
{
    friend Tree<NodeDataType, EdgeDataType>;

public:
    TreeNode() : link_(nullptr) {}

    // avoid copy constructor and assignment operator.
    // creating copies is maintained by Tree only.
    TreeNode (const TreeNode&) = delete;
    TreeNode& operator = (const TreeNode&) = delete;

    // -----------------------------------------------------
    //     Typedefs
    // -----------------------------------------------------

    typedef Tree    <NodeDataType, EdgeDataType> TreeType;
    typedef TreeLink<NodeDataType, EdgeDataType> LinkType;
    typedef TreeNode<NodeDataType, EdgeDataType> NodeType;
    typedef TreeEdge<NodeDataType, EdgeDataType> EdgeType;

    // -----------------------------------------------------
    //     Accessors
    // -----------------------------------------------------

    /**
     * @brief Returns the link of this node that points towards the root.
     */
    inline LinkType* PrimaryLink() const
    {
        return link_;
    }

    /**
     * @brief Returns the link of this node that points towards the root.
     *
     * This is just an alias for PrimaryLink(), that is shorter to use when needed
     * frequently in an algorithm.
     */
    inline LinkType* Link() const
    {
        return link_;
    }

    // -----------------------------------------------------
    //     Iterators
    // -----------------------------------------------------

    typedef TreeNodeIteratorLinks<      LinkType*,       NodeType*,       EdgeType*>      IteratorLinks;
    typedef TreeNodeIteratorLinks<const LinkType*, const NodeType*, const EdgeType*> ConstIteratorLinks;

    inline IteratorLinks BeginLinks()
    {
        return IteratorLinks(link_);
    }

    inline IteratorLinks EndLinks()
    {
        return IteratorLinks(nullptr);
    }

    inline ConstIteratorLinks BeginLinks() const
    {
        return ConstIteratorLinks(link_);
    }

    inline ConstIteratorLinks EndLinks() const
    {
        return ConstIteratorLinks(nullptr);
    }

    // -----------------------------------------------------
    //     Member Functions
    // -----------------------------------------------------

    /**
     * @brief Returns the index of this Link.
     */
    inline size_t Index() const
    {
        return index_;
    }

    int  Rank() const;
    bool IsLeaf() const;
    bool IsInner() const;

    //~ /** True if the node is the root, false otherwise. */
    //~ bool IsRoot();

    //~ /** Depth of the node in the tree, i.e. its distance from root. */
    //~ int Depth();

    std::string Dump() const;

    // -----------------------------------------------------
    //     Member Variables
    // -----------------------------------------------------

    //~ NodeDataType data;

// TODO !!! make protected again, and use some other mechanism for setting the members !!!
//~ protected:

    size_t index_;

    LinkType* link_;
};

// =============================================================================
//     Iterator Links
// =============================================================================

template <typename LinkPointerType, typename NodePointerType, typename EdgePointerType>
class TreeNodeIteratorLinks
{
public:
    // -----------------------------------------------------
    //     Typedefs
    // -----------------------------------------------------

    typedef TreeNodeIteratorLinks<LinkPointerType, NodePointerType, EdgePointerType> self_type;
    typedef std::forward_iterator_tag iterator_category;

    // -----------------------------------------------------
    //     Constructor
    // -----------------------------------------------------

    TreeNodeIteratorLinks (LinkPointerType link) :
    link_(link), start_(link)
    {}

    // -----------------------------------------------------
    //     Operators
    // -----------------------------------------------------

    inline self_type operator ++ ()
    {
        link_ = link_->Next();
        if (link_ == start_) {
            link_ = nullptr;
        }
        return *this;
    }

    inline self_type operator ++ (int)
    {
        self_type tmp = *this;
        ++(*this);
        return tmp;
    }

    inline bool operator == (const self_type &other) const
    {
        return other.link_ == link_;
    }

    inline bool operator != (const self_type &other) const
    {
        return !(other == *this);
    }

    // -----------------------------------------------------
    //     Members
    // -----------------------------------------------------

    inline LinkPointerType Link()
    {
        return link_;
    }

    inline NodePointerType Node()
    {
        return link_->Node();
    }

    inline EdgePointerType Edge()
    {
        return link_->Edge();
    }

    inline LinkPointerType StartLink()
    {
        return start_;
    }

protected:
    LinkPointerType link_;
    LinkPointerType start_;
};

} // namespace genesis

// =============================================================================
//     Inclusion of the implementation
// =============================================================================

// This is a class template, so do the inclusion here.
#include "tree/tree_node.tpp"

#endif // include guard
