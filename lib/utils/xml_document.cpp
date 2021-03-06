/**
 * @brief Implementation of XML document functions.
 *
 * @file
 * @ingroup utils
 */

#include "utils/xml_document.hpp"

#include "utils/logging.hpp"

namespace genesis {

// =============================================================================
//     Converter Functions
// =============================================================================

/**
 * @brief Converts a pointer to an XmlValue to XmlComment if appropriate.
 *
 * Triggers a warning and returns a nullptr if the dynamic type of the object is not actually
 * an XmlComment object.
 */
const XmlComment* XmlValueToComment (const XmlValue* v)
{
    if (v->type() != XmlValue::kComment) {
        LOG_WARN << "Invalid conversion from XmlValue::" << v->TypeToString()
                 << " to XmlValue::Comment.";
        return nullptr;
    }
    return static_cast<XmlComment*> (const_cast<XmlValue*> (v));
}

/**
 * @brief Converts a pointer to an XmlValue to XmlMarkup if appropriate.
 *
 * Triggers a warning and returns a nullptr if the dynamic type of the object is not actually
 * an XmlMarkup object.
 */
const XmlMarkup* XmlValueToMarkup (const XmlValue* v)
{
    if (v->type() != XmlValue::kMarkup) {
        LOG_WARN << "Invalid conversion from XmlValue::" << v->TypeToString()
                 << " to XmlValue::Markup.";
        return nullptr;
    }
    return static_cast<XmlMarkup*> (const_cast<XmlValue*> (v));
}

/**
 * @brief Converts a pointer to an XmlValue to XmlElement if appropriate.
 *
 * Triggers a warning and returns a nullptr if the dynamic type of the object is not actually
 * an XmlElement object.
 */
const XmlElement* XmlValueToElement (const XmlValue* v)
{
    if (v->type() != XmlValue::kElement) {
        LOG_WARN << "Invalid conversion from XmlValue::" << v->TypeToString()
                 << " to XmlValue::Element.";
        return nullptr;
    }
    return static_cast<XmlElement*> (const_cast<XmlValue*> (v));
}

} // namespace genesis
