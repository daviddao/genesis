/**
 * @brief Implementation of Lexer functions.
 *
 * @file
 * @ingroup utils
 */

#include "utils/lexer.hpp"

#include <assert.h>
#include <cstring>
#include <stack>
#include <string>

#include "utils/logging.hpp"
#include "utils/utils.hpp"

namespace genesis {

// =============================================================================
//     Process
// =============================================================================

/**
 * @brief Shortcut function that reads the contents of a file and then calls ProcessString().
 *
 * If the file does not exist, a warning is triggered and false returned. Otherwise, the result
 * of ProcessString() is returned.
 */
bool Lexer::ProcessFile(const std::string& fn)
{
    if (!FileExists(fn)) {
        LOG_WARN << "File '" << fn << "' does not exist.";
        return false;
    }
    return ProcessString(FileRead(fn));
}

/**
 * @brief Process a string and store the resulting tokens in this Lexer object.
 *
 * This function clears the token list stored for this object and fills it
 * with the results of processing the given string. This process analyzes and
 * splits the string into different tokens. For the types of tokens being
 * extracted, see LexerToken; for accessing the results, see Lexer.
 *
 * Returns true iff successful. In case an error is encountered while analyzing
 * the text, this functions returns false and the last token will be of type
 * LexerTokenType::kError, with the value being an error message describing the
 * type of error.
 *
 * Common usage:
 *
 *     Lexer l;
 *     l.ProcessString("tree(some:0.5,items:0.3);")
 *     if (l.empty()) {
 *         LOG_WARN << "Lexer is empty.";
 *     }
 *     if (l.HasError()) {
 *         LexerToken b = l.back();
 *         LOG_WARN << "Lexing error "
 *             << " at " << b.line() << ":" << b.column()
 *             << " with message " << b.value();
 *     }
 *     ... process the tokens ...
 *
 * The additional option `stepwise` will not scan the entire string, but only the
 * first element of it. In order to get more tokens manually, ProcessStep() has to be called until
 * it returns false. However, this option is usually used in combination with a producing
 * LexerIterator, that will automatically call ProcessStep() whenever tokens are needed.
 * See JsonProcessor::FromString() for an example.
 */
bool Lexer::ProcessString(const std::string& text, bool stepwise)
{
    Init(text);

    // if we want stepwise lexing, just do the first step.
    if (stepwise) {
        return ProcessStep();
    }

    // if not, do steps till the end.
    while (!IsEnd()) {
        if (!ProcessStep()) {
            return tokens_.empty();
        }
    }

    return true;
}

/**
 * @brief Processes one step of the lexing.
 *
 * This might produce more than one token, as comments and whitespaces are treaded specially.
 *
 * As stated in the description of this Lexer class, the class is meant to be
 * derived for concrete types of lexers. Thus, here are some comments about the
 * working of this function:
 *
 * For most types of text files, the first character of each token determines
 * the type of the token (for example, a digit almost always leads to a number
 * token). This is why we use a list telling us which char leads to which token
 * type. This list is a speedup, because using it, we do not need to try every
 * scanner (for numbers, symbols, strings, etc) at the beginning of each new
 * token, but simply do a lookup to find out "this char means we have to use
 * this scanner now". (Also see GetCharType and SetCharType for this.)
 *
 * This does not mean that any char of a given type can only appear in
 * tokens of that type. For example, typically a symbol can start with
 * letters, but then contain numbers later, too. Thus, although a digit
 * is of type kNumber, it can end up in a symbol token, depending on the
 * context.
 *
 * Instead, the char type here is meant as a speedup for finding the right
 * scanner when lexing the text: whenever one scanner finishes, the next
 * char is inspected and depending on its type, a fitting scanner is
 * activated (for digits the number scanner and so on).
 *
 * This technique will not work if finding the correct scanner depends on
 * more than the first character of the token. For example, comments usually
 * start with a more complex sequence ("//" or even "<!--"), which is why
 * they are specially treaded in the ProcessString function.
 *
 * So, in situations, where the type of the next token cannot be determined from
 * its first character (except comments), ProcessString has to be overridden in the
 * derived class in order to do some other checking methods to determine the
 * correct scanner. In the new ProcessString function, first call Init to reset all
 * internal variables. Also see ScanUnknown for some important information.
 */
bool Lexer::ProcessStep()
{
    if (IsEnd()) {
        return false;
    }

    // scan arbitrary amount of interleaved whitespace and comments
    while (ScanWhitespace() || ScanComment()) continue;

    // check if whitespace or comment scanner yielded an error
    if (!tokens_.empty() && tokens_.back().IsError()) {
        return false;
    }

    // check if still not at end after whitespaces and comments
    if (IsEnd()) {
        return false;
    }

    // check if current char is an error char
    LexerTokenType t = GetCharType();
    if (t == LexerTokenType::kError) {
        PushToken(LexerTokenType::kError, GetPosition(), "Invalid character.");
        return false;
    }

    // start the actual scanners depending on the first char
    switch (t) {
        case LexerTokenType::kSymbol:
            ScanSymbol();
            break;
        case LexerTokenType::kNumber:
            ScanNumber();
            break;
        case LexerTokenType::kString:
            ScanString();
            break;
        case LexerTokenType::kBracket:
            ScanBracket();
            break;
        case LexerTokenType::kOperator:
            ScanOperator();
            break;
        case LexerTokenType::kTag:
            ScanTag();
            break;
        case LexerTokenType::kUnknown:
            ScanUnknown();
            break;

        case LexerTokenType::kWhite:
        case LexerTokenType::kComment:
            // this case happens if ScanWhitespace or ScanComment are wrongly
            // implemented in a derived class, so that they return false
            // although they should scan something or create an error token.
            // this is particularly the case when ScanComment does not create
            // an error token when a close-comment-char appears without a
            // open-comment-char first, e.g. "]".
            assert(false);
        case LexerTokenType::kError:
            // this cannot happen, as we already covered the case before the switch
            assert(false);
        default:
            // this case only happens if someone added a new element to the enum.
            // make sure that we do not forget to include it in this switch!
            assert(false);
    }

    // check if the scanners produced an error
    if (tokens_.empty() || tokens_.back().IsError()) {
        return false;
    }
    return true;
}

//~ void Lexer::StripComments ()
//~ {
//~
//~ }

// =============================================================================
//     Scanners
// =============================================================================

/**
 * @brief Scans a range between two strings.
 *
 * If the current position in the text starts with the value of `from`, this
 * scanner continues in the text until the value of `to` is found (or the end of
 * the text). In case of success (both `from` and `to` were found), it returns true,
 * false otherwise.
 *
 * If successfull, the internal position of the Lexer is moved to the char after
 * the `to` string. Thus, the start position has to be stored before calling this function.
 *
 * Example:
 *
 *     size_t start = GetPosition();
 *     bool   found = ScanFromTo("{", "}");
 *     if (found) {
 *         PushToken(LexerTokenType::kTag, start+1, GetPosition()-1);
 *     }
 *
 * This scans between two curly brackets and if found, produces a token (that does not contain
 * the brackets themselves).
 */
bool Lexer::ScanFromTo (const char* from, const char* to)
{
    // first check if the current position actually contains the "from" string
    if (IsEnd() || strncmp(from, text_+itr_, strlen(from)) != 0) {
        return false;
    }

    // if so, move as many chars forward. we have to split this from the
    // checking, because we do not want to change itr_ in case it is not a
    // match. also, calling NextChar here ensures integrity of the line
    // counting.
    for (size_t i = 0; i < strlen(from); ++i) {
        NextChar();
    }

    // now try to find the "to" string
    while (!IsEnd() && strncmp(to, text_+itr_, strlen(to)) != 0) {
        NextChar();
    }

    // if the "to" string was not found before the end of the text, we are done
    if (IsEnd()) {
        return false;
    }

    // "to" string was found. move as many chars forward.
    for (size_t i = 0; i < strlen(to); ++i) {
        NextChar();
    }
    return true;
}

/**
 * @brief Scans the text as long as the current char is of type kUnknown.
 *
 * It is possible that this function has to be overridden in case that ProcessString()
 * is overridden as well. This is because of the following:
 *
 * Overriding ProcessString means that the approach of determining the correct scanner
 * for a token from its first character does not work. Thus, all characters that
 * do not indicate a scanner will usually be set to unknown in order to treat
 * them specially in ProcessString. This means that scanning unknown chars should stop
 * after each char, return to ProcessString(), and check again if now it is time to
 * activate a special scanner. Thus, in this case, ScanUnknown should only scan
 * one character at a time.
 */
inline bool Lexer::ScanUnknown()
{
    size_t start = GetPosition();
    while (!IsEnd() && GetCharType() == LexerTokenType::kUnknown) {
        NextChar();
    }
    PushToken(LexerTokenType::kUnknown, start, GetPosition());
    return true;
}

/*
// version that only scans one char at a time. useful for allowing e.g.
// complex comment markers, that would otherwise be turned into unknown tokens.
inline bool Lexer::ScanUnknown()
{
    PushToken(LexerTokenType::kUnknown, GetPosition(), GetPosition() + 1);
    NextChar();
    return true;
}
*/

/**
 * @brief Scan for whitespace.
 *
 * Returns true iff whitespace was found.
 */
inline bool Lexer::ScanWhitespace()
{
    bool   found = false;
    size_t start = GetPosition();

    while (!IsEnd() && GetCharType() == LexerTokenType::kWhite) {
        NextChar();
        found = true;
    }
    if (include_whitespace && found) {
        PushToken(LexerTokenType::kWhite, start, GetPosition());
    }
    return found;
}

/**
 * @brief Scans for comments in different formats.
 *
 * In the base class, this functions simply returns false. In order to scan for
 * actual comments, it has to be overridden. A typical implementation can be
 * found in NewickLexer::ScanComment().
 */
inline bool Lexer::ScanComment()
{
    return false;
}

/**
 * @brief Scan a symbol.
 *
 * In the base class, the functions scans as long as the LexerToken type of the
 * current char is kSymbol. If it is necessary to continue scanning for other
 * types of chars once a symbol has started, this function has to be overridden.
 * This is particularly the case for digits: typically, symbols cannot start
 * with a digit, but once a symbol has started, digits can occur.
 *
 * Returns true, as symbols cannot be malformatted.
 */
inline bool Lexer::ScanSymbol()
{
    size_t start = GetPosition();
    while (!IsEnd() && GetCharType() == LexerTokenType::kSymbol) {
        NextChar();
    }
    PushToken(LexerTokenType::kSymbol, start, GetPosition());
    return true;
}

/**
 * @brief Scans a number.
 *
 * Following formats are supported:
 *
 *     [+-]123[.456][eE[+-]789]
 *
 * Returns true iff valid number was found.
 */
inline bool Lexer::ScanNumber()
{
    size_t start   = GetPosition();
    bool   found_d = false; // found a dot
    bool   found_e = false; // found the letter e

    // encountered an error while scanning. this happens if we reach the end
    // of the number before seeing any number content (digit, dot, etc)
    bool   err     = false;

    // scan
    while(!IsEnd()) {
        if(CharIsDigit(GetChar())) {
            // nothing to do
        } else if (GetChar() == '.') {
            // do not allow more than one dot, require a number after the dot
            // (if not, treat it as the end of the number, stop scanning)
            if (
                found_d
                || IsEnd(+1) || !CharIsDigit(GetChar(+1))
            ) {
                break;
            }
            found_d = true;
        } else if (CharMatch(GetChar(), 'e')) {
            // do not allow more than one e (treat the second one as the end of
            // the number and stop scannning).
            // also, require a number or sign before and after the first e. if not,
            // treat it also as the end of the number and stop scanning
            if (
                found_e
                || GetPosition() == 0 || !CharIsDigit(GetChar(-1))
                || IsEnd(+1)
                || (!CharIsDigit(GetChar(+1)) && !CharIsSign(GetChar(+1)))
            ) {
                err = (GetPosition() == start);
                break;
            }
            found_e = true;
        } else if (CharIsSign(GetChar())) {
            // conditions for when a sign is valid:
            //   - it is at the beginning of the token and followed by digits
            //   - it comes immediately after the e and is followed by digits
            // --> when neither is fullfilled, treat it as the end of the
            //     number, stop scanning
            if (
                !(
                    GetPosition() == start
                    && !IsEnd(+1) && CharIsDigit(GetChar(+1))
                ) &&
                !(
                    found_e && CharMatch(GetChar(-1), 'e')
                    && !IsEnd(+1) && CharIsDigit(GetChar(+1))
                )
            ) {
                err = (GetPosition() == start);
                break;
            }
        } else {
            err = (GetPosition() == start);
            break;
        }
        NextChar();
    }

    // create result
    if (err) {
        PushToken(LexerTokenType::kError, GetPosition(), "Malformed number.");
        return false;
    } else {
        PushToken(LexerTokenType::kNumber, start, GetPosition());
        return true;
    }
}

/**
 * @brief Scan a string.
 *
 * A string is usually enclosed either in 'abc' or in "def" (this depends on
 * which chars have been set to LexerTokenType kString in the derived class).
 * Within a string, any character is allowed. See use_string_escape and
 * use_string_doubled_quotes for options modifying the behaviour of this
 * function.
 *
 * Returns true iff the string is finished with the correct quotation mark.
 */
inline bool Lexer::ScanString()
{
    // skip the first quotation mark, save its value for later comparision
    // so that the string ends with the same type of mark
    char qmark = GetChar();
    NextChar();
    if (IsEnd()) {
        PushToken(LexerTokenType::kError, GetPosition()-1, "Malformed string.");
        return false;
    }

    size_t start = GetPosition();
    bool jump    = false; // marks when we jumped over an escape or quote mark
    bool found_e = false; // found an escape sequence
    bool found_q = false; // found a doubled qutation mark ""

    // scan
    while (!IsEnd()) {
        // if we find a backslash and use escape characters, we skip the
        // backslash and the following char. they will then be de-escaped after
        // the end of the string is reached.
        if (GetChar() == '\\' && use_string_escape) {
            jump    = true;
            found_e = true;
            NextChar();
            NextChar();
            continue;
        }
        if (
            GetChar() == qmark && GetChar(+1) == qmark
            && use_string_doubled_quotes
        ) {
            jump    = true;
            found_q = true;
            NextChar();
            NextChar();
            continue;
        }
        // check if we reached the end of the string
        if (GetChar() == qmark) {
            NextChar();
            break;
        }
        jump = false;
        NextChar();
    }

    // reached end of text before ending quotation mark.
    // we need to check jump here, because an escape sequence or doubled quote mark
    // (those are the situations when jump=true) can look like the string ending,
    // but actually isn't. this means, if jump is true, we left the loop right after such a case
    // (escape sequence or doubled qmark), and are now at the end of the text, without having
    // seen the end qmark. e.g.: "hello world\"
    if (jump || (IsEnd() && !(GetChar(-1) == qmark))) {
        PushToken(LexerTokenType::kError, start-1, "Malformed string.");
        return false;
    }

    // de-escape the string (transform backslash-escaped chars)
    std::string res = GetSubstr(start, GetPosition()-1);
    if (found_e && use_string_escape) {
        res = StringDeescape(res);
    }

    // transform doubled qmarks like "" into single ones like "
    if (found_q && use_string_doubled_quotes) {
        res = StringReplaceAll(
            res,
            std::string(2, qmark),
            std::string(1, qmark)
        );
    }

    // if needed, add qmarks again
    if (!trim_quotation_marks) {
        res = qmark + res + qmark;
    }

    // create result
    PushToken(LexerTokenType::kString, start-1, res);
    return true;
}

/**
 * @brief Scans a single operators.
 *
 * If the operator is a sign and the next char in the text is a digit,
 * and glue_sign_to_number is set, we scan it as a number. For example, the
 * sequence
 *
 *     a+=-3;
 *
 * will result in four tokens: Symbol "a", Operators "+" and "=" and Number "-3",
 * instead of five tokens with three consequtive operators.
 *
 * This allows to use signs (+-) as operator chars in cases where they are not
 * followed by a number.
 */
inline bool Lexer::ScanOperator()
{
    // if the operator is a sign followed by a number, scan it as a number
    if (
        CharIsSign(GetChar()) && glue_sign_to_number &&
        !IsEnd(+1) && CharIsDigit(GetChar(+1))
    ) {
        return ScanNumber();
    }

    PushToken(LexerTokenType::kOperator, GetPosition(), GetPosition()+1);
    NextChar();
    return true;
}

/**
 * @brief Scans a single bracket.
 *
 * Returns true.
 */
inline bool Lexer::ScanBracket()
{
    PushToken(LexerTokenType::kBracket, GetPosition(), GetPosition()+1);
    NextChar();
    return true;
}

/**
 * @brief Scans a tag.
 *
 * Tags can be anything from a single char to a complete string that has to
 * be processed further. It might be neccessary to run a second lexer step using
 * another lexer class on the inner part of the tag, for example when scanning
 * xml text:
 *
 *     <tag key1="x" key2-="y">
 *
 * may result as one tag (with or without the enclosing <>, depending on the
 * specific implementation and needs). Then, a second lexer can decompose this into
 *
 *     tag key1 = x key2 = y
 *
 * In the base class, this functions simply returns false. In order to scan for
 * actual comments, it has to be overridden. A typical implementation can be
 * found in NewickLexer::ScanComment().
 */
inline bool Lexer::ScanTag()
{
    return false;
}

// =============================================================================
//     Iterator
// =============================================================================

/**
 * @brief Returns an iterator to the beginning of the token list.
 *
 * This is used for the iterator and also allows to use range based
 * looping over the tokens:
 *
 *     Lexer l;
 *     for (LexerToken& t : l) {
 *         std::cout << t.value() << std::endl;
 *     }
 * %
 */
inline Lexer::iterator Lexer::begin()
{
    return iterator(*this, 0);
}

/**
 * @brief Returns an iterator to the end of the token list.
 */
inline Lexer::iterator Lexer::end()
{
    return iterator(*this, -1);
}

// =============================================================================
//     Helper functions
// =============================================================================

/**
 * @brief Checkes whether the bracket tokes are validly nested.
 *
 * In order to be valid, every opening bracket must be matched with a
 * corresponding closing bracket, and their order has to be correct.
 *
 * This function checks the following types of brackets:
 *
 *     () [] {} <>
 *
 * It is not particularly useful for xml, as there it is also important to use
 * closing tags like `<xml> ... </xml>`.
 */
bool Lexer::ValidateBrackets() const
{
    std::stack<char> stk;
    for (LexerToken t : tokens_) {
        if (!t.IsBracket()) {
            continue;
        }

        char c = t.value()[0];
        if (c == '(') stk.push(')');
        if (c == '[') stk.push(']');
        if (c == '{') stk.push('}');
        if (c == '<') stk.push('>');

        if (c == ')' || c == ']' || c == '}' || c == '>') {
            if (stk.empty() || c != stk.top()) {
                return false;
            } else {
                stk.pop();
            }
        }
    }
    return true;
}

/**
 * @brief Returns a listing of the parse result in readable form.
 */
std::string Lexer::Dump() const
{
    std::string res;
    for (size_t i = 0; i < tokens_.size(); i++) {
        LexerToken t = tokens_[i];
        char out[30];
        sprintf(out, "[%03d] @%03d:%03d %10s : ",
            static_cast<unsigned int>(i),
            static_cast<unsigned int>(t.line()),
            static_cast<unsigned int>(t.column()),
            t.TypeToString().c_str()
        );
        res += out + t.value() + '\n';
    }
    return res;
}

} // namespace genesis
