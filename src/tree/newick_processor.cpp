/**
 * @brief Implementation of functions for reading and writing Newick files.
 *
 * For reasons of readability, in this implementation file, the template data types
 * NodeDataType and EdgeDataType are abbreviated using NDT and EDT, respectively.
 *
 * @file
 * @ingroup tree
 */

#include "tree/newick_processor.hpp"

namespace genesis {

// =============================================================================
//     Printing
// =============================================================================

bool NewickProcessor::print_names          = true;
bool NewickProcessor::print_branch_lengths = false;
bool NewickProcessor::print_comments       = false;
bool NewickProcessor::print_tags           = false;

// TODO this is a quick and dirty (=slow) solution...
std::string NewickProcessor::ToStringRec(NewickBroker& broker, size_t pos)
{
    // check if it is a leaf, stop recursion if so.
    if (broker[pos]->rank() == 0) {
        return ElementToString(broker[pos]);
    }

    // recurse over all children of the current node. while doing so, build a stack of the resulting
    // substrings in reverse order. this is because newick stores the nodes kind of "backwards",
    // by starting at a leaf node instead of the root.
    std::deque<std::string> children;
    for (size_t i = pos + 1; i < broker.size() && broker[i]->depth > broker[pos]->depth; ++i) {
        // skip if not immediate children (those will be called in later recursion steps)
        if (broker[i]->depth > broker[pos]->depth + 1) {
            continue;
        }

        // do the recursion step for this child, add the result to a stack
        children.push_front(ToStringRec(broker, i));
    }

    // build the string by iterating the stack
    std::ostringstream out;
    out << "(";
    for (size_t i = 0; i < children.size(); ++i) {
        if (i>0) {
            out << ",";
        }
        out << children[i];
    }
    out << ")" << ElementToString(broker[pos]);
    return out.str();
}

std::string NewickProcessor::ElementToString(NewickBrokerElement* bn)
{
    std::string res = "";
    if (print_names) {
        res += StringReplaceAll(bn->name, " ", "_");
    }
    if (print_branch_lengths) {
        res += ":" + std::to_string(bn->branch_length);
    }
    if (print_comments) {
        for (std::string c : bn->comments) {
            res += "[" + c + "]";
        }
    }
    if (print_tags) {
        for (std::string t : bn->tags) {
            res += "{" + t + "}";
        }
    }
    return res;
}

} // namespace genesis
