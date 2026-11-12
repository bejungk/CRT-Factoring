#include "xlsx_export.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace factoring {
namespace {

using Bytes = std::string;

/**
 * @brief Append a 16-bit integer to a byte string in little-endian order.
 * @param out Destination byte string.
 * @param value Integer to append.
 */
void append_u16(Bytes& out, std::uint16_t value)
{
    // ZIP format integers are encoded least-significant byte first.
    out.push_back(static_cast<char>(value & 0xff));
    out.push_back(static_cast<char>((value >> 8) & 0xff));
}

/**
 * @brief Append a 32-bit integer to a byte string in little-endian order.
 * @param out Destination byte string.
 * @param value Integer to append.
 */
void append_u32(Bytes& out, std::uint32_t value)
{
    for (unsigned int shift = 0; shift < 32; shift += 8) {
        out.push_back(static_cast<char>((value >> shift) & 0xff));
    }
}

/**
 * @brief Compute the ZIP CRC-32 checksum for a byte sequence.
 * @param data Bytes to checksum.
 * @return CRC-32 checksum.
 */
std::uint32_t crc32(std::string_view data)
{
    // Compute CRC-32 (IEEE polynomial) for an uncompressed ZIP member.
    std::uint32_t crc = 0xffffffff;
    for (unsigned char byte : data) {
        crc ^= byte;
        for (int bit = 0; bit < 8; ++bit) {
            crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
        }
    }
    return ~crc;
}

/**
 * @brief Escape XML metacharacters in text content.
 * @param value Unescaped text.
 * @return XML-safe text.
 */
std::string xml_escape(std::string_view value)
{
    // Escape XML markup characters in spreadsheet cell text.
    std::string escaped;
    for (char ch : value) {
        switch (ch) {
        case '&': escaped += "&amp;"; break;
        case '<': escaped += "&lt;"; break;
        case '>': escaped += "&gt;"; break;
        case '"': escaped += "&quot;"; break;
        case '\'': escaped += "&apos;"; break;
        default: escaped += ch; break;
        }
    }
    return escaped;
}

/**
 * @brief Convert a one-based column index to Excel letter notation.
 * @param index One-based column index.
 * @return Column letters such as A, Z, or AA.
 */
std::string column_name(std::size_t index)
{
    // Convert a one-based column index to Excel's A..Z, AA.. naming scheme.
    std::string name;
    while (index != 0) {
        --index;
        name.insert(name.begin(), static_cast<char>('A' + index % 26));
        index /= 26;
    }
    return name;
}

using Sheet = std::vector<std::vector<std::string>>;

/**
 * @brief Serialize a table as a minimal OOXML worksheet.
 * @param rows Rows of cell text.
 * @return Worksheet XML document.
 */
std::string worksheet_xml(const Sheet& rows)
{
    // Emit a minimal worksheet whose cells all hold inline strings.
    std::string xml = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
                      "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\"><sheetData>";
    for (std::size_t r = 0; r < rows.size(); ++r) {
        xml += "<row r=\"" + std::to_string(r + 1) + "\">";
        for (std::size_t c = 0; c < rows[r].size(); ++c) {
            const std::string ref = column_name(c + 1) + std::to_string(r + 1);
            xml += "<c r=\"" + ref + "\" t=\"inlineStr\"><is><t xml:space=\"preserve\">" +
                   xml_escape(rows[r][c]) + "</t></is></c>";
        }
        xml += "</row>";
    }
    xml += "</sheetData></worksheet>";
    return xml;
}

struct ZipEntry {
    std::string name;
    std::string data;
};

/**
 * @brief Write ZIP entries as an uncompressed archive to disk.
 * @param path Destination archive path.
 * @param entries Named archive members and their contents.
 */
void write_zip(const std::string& path, const std::vector<ZipEntry>& entries)
{
    // Package the OOXML parts in an uncompressed ZIP archive, including its
    // local headers, central directory, and end-of-directory record.
    Bytes archive;
    struct CentralEntry {
        const ZipEntry* entry;
        std::uint32_t crc;
        std::uint32_t offset;
    };
    std::vector<CentralEntry> central;
    for (const ZipEntry& entry : entries) {
        const std::uint32_t crc = crc32(entry.data);
        const std::uint32_t offset = static_cast<std::uint32_t>(archive.size());
        append_u32(archive, 0x04034b50);
        append_u16(archive, 20);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u32(archive, crc);
        append_u32(archive, static_cast<std::uint32_t>(entry.data.size()));
        append_u32(archive, static_cast<std::uint32_t>(entry.data.size()));
        append_u16(archive, static_cast<std::uint16_t>(entry.name.size()));
        append_u16(archive, 0);
        archive += entry.name;
        archive += entry.data;
        central.push_back({&entry, crc, offset});
    }

    const std::uint32_t central_offset = static_cast<std::uint32_t>(archive.size());
    for (const CentralEntry& item : central) {
        append_u32(archive, 0x02014b50);
        append_u16(archive, 20);
        append_u16(archive, 20);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u32(archive, item.crc);
        const auto size = static_cast<std::uint32_t>(item.entry->data.size());
        append_u32(archive, size);
        append_u32(archive, size);
        append_u16(archive, static_cast<std::uint16_t>(item.entry->name.size()));
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u16(archive, 0);
        append_u32(archive, 0);
        append_u32(archive, item.offset);
        archive += item.entry->name;
    }
    const std::uint32_t central_size =
        static_cast<std::uint32_t>(archive.size()) - central_offset;
    append_u32(archive, 0x06054b50);
    append_u16(archive, 0);
    append_u16(archive, 0);
    append_u16(archive, static_cast<std::uint16_t>(entries.size()));
    append_u16(archive, static_cast<std::uint16_t>(entries.size()));
    append_u32(archive, central_size);
    append_u32(archive, central_offset);
    append_u16(archive, 0);

    std::ofstream file(path, std::ios::binary);
    if (!file || !file.write(archive.data(), static_cast<std::streamsize>(archive.size()))) {
        throw std::runtime_error("Could not write XLSX file: " + path);
    }
}

} // namespace

/**
 * @brief Export summary, factor-base, and relation tables to an XLSX file.
 * @param path Destination file path.
 * @param number Integer N represented by the workbook.
 * @param bound Reported factor-base bound value.
 * @param factor_base Selected primes and modular roots.
 * @param relations CRT relations to include.
 */
void export_results_xlsx(
    const std::string& path, const mpz_class& number, std::uint32_t bound,
    const std::vector<factor_base_entry>& factor_base,
    const std::vector<crt_relation>& relations)
{
    // Construct three worksheet tables and the package metadata linking them.
    Sheet summary{{"N", number.get_str()}, {"B (floor(log2(N)))", std::to_string(bound)},
                  {"Relation count", std::to_string(relations.size())}};
    Sheet base{{"Prime", "Root 1", "Root 2"}};
    for (const factor_base_entry& entry : factor_base) {
        base.push_back({entry.prime.get_str(), entry.root1.get_str(), entry.root2.get_str()});
    }
    Sheet relation_rows{{"Combination", "x", "Multiplier g", "Modulus",
                         "Prime", "Selected root",
                         "Exponent", "Remaining cofactor", "Cofactor is prime"}};
    for (std::size_t r = 0; r < relations.size(); ++r) {
        const crt_relation& relation = relations[r];
        const bool cofactor_is_prime =
            relation.remaining_cofactor > 1 &&
            mpz_probab_prime_p(relation.remaining_cofactor.get_mpz_t(), 25) != 0;
        for (std::size_t i = 0; i < factor_base.size(); ++i) {
            relation_rows.push_back({std::to_string(r + 1), relation.x.get_str(),
                                     relation.multiplier.get_str(),
                                     relation.modulus.get_str(),
                                     factor_base[i].prime.get_str(),
                                     relation.prime_roots[i].second.get_str(),
                                     std::to_string(relation.exponent_vector[i]),
                                     relation.remaining_cofactor.get_str(),
                                     cofactor_is_prime ? "Yes" : "No"});
        }
    }

    const std::string content_types =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet2.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet3.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "</Types>";
    const std::string root_rels =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
        "</Relationships>";
    const std::string workbook =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><sheets>"
        "<sheet name=\"Summary\" sheetId=\"1\" r:id=\"rId1\"/>"
        "<sheet name=\"Factor Base\" sheetId=\"2\" r:id=\"rId2\"/>"
        "<sheet name=\"Relations\" sheetId=\"3\" r:id=\"rId3\"/>"
        "</sheets></workbook>";
    const std::string workbook_rels =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
        "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet2.xml\"/>"
        "<Relationship Id=\"rId3\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet3.xml\"/>"
        "</Relationships>";

    write_zip(path, {{"[Content_Types].xml", content_types},
                     {"_rels/.rels", root_rels},
                     {"xl/workbook.xml", workbook},
                     {"xl/_rels/workbook.xml.rels", workbook_rels},
                     {"xl/worksheets/sheet1.xml", worksheet_xml(summary)},
                     {"xl/worksheets/sheet2.xml", worksheet_xml(base)},
                     {"xl/worksheets/sheet3.xml", worksheet_xml(relation_rows)}});
}

/**
 * @brief Export collected cofactor experiment results to an XLSX file.
 * @param path Destination XLSX file.
 * @param results Cofactor observations from the experiment.
 */
void export_cofactors_xlsx(
    const std::string& path,
    const std::vector<cofactor_experiment_result>& results)
{
    Sheet rows{{"Factor-base stage", "Combination", "Prime A", "Prime B", "N",
                "M", "x", "x bits", "Cofactor", "Cofactor bits",
                "Cofactor primality", "Duplicate cofactor?",
                "Duplicate relation?",
                "Known prime factors (prime^exponent)"}};
    for (const cofactor_experiment_result& result : results) {
        const int primality_result =
            mpz_probab_prime_p(result.cofactor.get_mpz_t(), 25);
        const char* primality_label =
            result.cofactor == 1 ? "Unit (1)"
            : primality_result == 0 ? "Composite"
            : primality_result == 2 ? "Prime"
                                    : "Probable prime";
        std::ostringstream known_factors;
        bool has_known_factor = false;
        const std::size_t exponent_count = std::min(
            result.exponent_primes.size(), result.exponent_vector.size());
        for (std::size_t i = 0; i < exponent_count; ++i) {
            if (result.exponent_vector[i] != 0) {
                if (has_known_factor) {
                    known_factors << " * ";
                }
                known_factors << result.exponent_primes[i].get_str() << '^'
                              << result.exponent_vector[i];
                has_known_factor = true;
            }
        }
        if (!has_known_factor) {
            known_factors << '1';
        }
        rows.push_back({std::to_string(result.experiment),
                        std::to_string(result.combination),
                        result.first_prime.get_str(), result.second_prime.get_str(),
                        result.number.get_str(), result.modulus.get_str(),
                        result.x.get_str(),
                        std::to_string(mpz_sizeinbase(result.x.get_mpz_t(), 2)),
                        result.cofactor.get_str(),
                        std::to_string(mpz_sizeinbase(result.cofactor.get_mpz_t(), 2)),
                        primality_label,
                        result.duplicate_cofactor ? "Yes" : "No",
                        result.duplicate_relation ? "Yes" : "No",
                        known_factors.str()});
    }

    const std::string content_types =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "</Types>";
    const std::string root_rels =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
        "</Relationships>";
    const std::string workbook =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\"><sheets>"
        "<sheet name=\"Cofactors\" sheetId=\"1\" r:id=\"rId1\"/>"
        "</sheets></workbook>";
    const std::string workbook_rels =
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
        "</Relationships>";

    write_zip(path, {{"[Content_Types].xml", content_types},
                     {"_rels/.rels", root_rels},
                     {"xl/workbook.xml", workbook},
                     {"xl/_rels/workbook.xml.rels", workbook_rels},
                     {"xl/worksheets/sheet1.xml", worksheet_xml(rows)}});
}

} // namespace factoring
