#ifndef GENESIS_TREE_NEWICKPROCESSOR_H_
#define GENESIS_TREE_NEWICKPROCESSOR_H_

/**
 * @brief
 *
 * @file
 * @ingroup tree
 */

#include <assert.h>
#include <string>

#include "utils/lexer.hpp"

namespace genesis {

// =============================================================================
//     Forward declarations
// =============================================================================

template <class NodeDataType, class EdgeDataType>
class  Tree;

class  NewickBroker;
struct NewickBrokerElement;

// =============================================================================
//     Newick Lexer
// =============================================================================

class NewickLexer : public Lexer
{
public:
    NewickLexer() {
        // set the special chars for newick trees
        SetCharType (LexerTokenType::kComment,  "[]");
        SetCharType (LexerTokenType::kTag,      "{}");
        SetCharType (LexerTokenType::kBracket,  "()");
        SetCharType (LexerTokenType::kOperator, ",;");

        // we use symbols and strings the same way here: both are labels for nodes, the first begin
        // called unquoted_label, the second quoted_label.
        SetCharType (LexerTokenType::kString,   "'");

        // the only numbers in newick are branch lengths, which are always introduced by a leading
        // colon, so we need only this here as starter for a number.
        SetCharType (LexerTokenType::kNumber,   ":");

        // this also allows (in accordance to the newick standard) to start a label with a digit.
        SetCharType (LexerTokenType::kSymbol,   "0123456789");

        // furthermore, set all remaining graphic chars to symbol so that they can be in a label.
        SetCharType (LexerTokenType::kSymbol,   "!\"#$%&*+-./<=>?@\\^_`|~");

        // set the flags as needed
        include_whitespace        = false;
        include_comments          = true;
        glue_sign_to_number       = false;
        trim_quotation_marks      = true;
        use_string_escape         = false;
        use_string_doubled_quotes = true;
    }

    /** @brief Determines whether to include tags in curly brackets. */
    bool include_tags = true;

protected:
    inline bool ScanComment()
    {
        if (GetChar() == ']') {
            PushToken(LexerTokenType::kError, GetPosition(), "Closing comment without opening it.");
            return false;
        }
        size_t start = GetPosition();
        bool   found = ScanFromTo("[", "]");
        if (!found && GetChar() == '[') {
            PushToken(LexerTokenType::kError, GetPosition(), "Comment not closed.");
            return false;
        }
        if (found && include_comments) {
            PushToken(LexerTokenType::kComment, start+1, GetPosition()-1);
        }
        return found;
    }

    inline bool ScanNumber()
    {
        // colon is the only char that is set as the beginning char for a number, so when we enter
        // this function, it must be the current char in the text
        assert(GetChar() == ':');

        // skip it and continue to scan as a normal number
        NextChar();
        return Lexer::ScanNumber();
    }

    inline bool ScanTag()
    {
        if (GetChar() == '}') {
            PushToken(LexerTokenType::kError, GetPosition(), "Closing tag without opening tag.");
            return false;
        }

        // curly brackets are the only chars that are used as tag in this lexer. we already checked
        // that the current char is not a closing bracket, so it must be an opening one.
        assert(GetChar() == '{');

        size_t start = GetPosition();
        bool   found = ScanFromTo("{", "}");
        if (!found) {
            PushToken(LexerTokenType::kError, start, "Opening tag without closing tag.");
            return false;
        }
        if (include_tags) {
            PushToken(LexerTokenType::kTag, start+1, GetPosition()-1);
        }
        return true;
    }
};

// =============================================================================
//     Newick Processor
// =============================================================================

class NewickProcessor
{
public:

    static std::string default_leaf_name;
    static std::string default_internal_name;
    static std::string default_root_name;

    static bool        use_default_names;

    // ---------------------------------------------------------------------
    //     Parsing
    // ---------------------------------------------------------------------

    template <class NodeDataType, class EdgeDataType>
    static bool FromFile   (const std::string fn,     Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static bool FromString (const std::string ts,     Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static bool FromLexer  (const NewickLexer& lexer, Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static void FromBroker (NewickBroker& broker,     Tree<NodeDataType, EdgeDataType>& tree);

    // ---------------------------------------------------------------------
    //     Printing
    // ---------------------------------------------------------------------

    static bool print_names;
    static bool print_branch_lengths;
    static bool print_comments;
    static bool print_tags;

    static int  precision;

    template <class NodeDataType, class EdgeDataType>
    static bool ToFile   (const std::string fn, const Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static void ToString (std::string& ts,      const Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static std::string ToString (               const Tree<NodeDataType, EdgeDataType>& tree);

    template <class NodeDataType, class EdgeDataType>
    static void ToBroker (NewickBroker& broker, const Tree<NodeDataType, EdgeDataType>& tree);

protected:
    static std::string ToStringRec(const NewickBroker& broker, size_t position);
    static std::string ElementToString(const NewickBrokerElement* bn);
};

} // namespace genesis

// =============================================================================
//     Inclusion of the implementation
// =============================================================================

// This class contains function templates, so do the inclusion here.
#include "tree/newick_processor.tpp"

#endif // include guard
