// Build with `make test`.
#include <algorithm>
#include <iostream>
#include "smith_waterman.h"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; std::cerr << "FAIL line " << __LINE__ << ": " #cond "\n"; } } while (0)

static std::string strip(std::string s) { s.erase(std::remove(s.begin(), s.end(), '-'), s.end()); return s; }

int main() {
    // Wikipedia example: match +3, mismatch -3, gap -2
    Scoring wiki{3, -3, -2};
    Alignment w = smithWaterman("GGTTGACTA", "TGTTACGG", wiki);
    CHECK(w.score == 13);
    CHECK(w.alignedA == "GTTGAC");
    CHECK(w.alignedB == "GTT-AC");

    // identical sequences align fully
    Alignment same = smithWaterman("ACGTACGT", "ACGTACGT");
    CHECK(same.score == 16);
    CHECK(same.alignedA == "ACGTACGT" && same.alignedB == "ACGTACGT");

    // no similarity -> empty alignment
    Alignment none = smithWaterman("AAAA", "TTTT");
    CHECK(none.score == 0 && none.alignedA.empty());

    // local alignment finds the embedded segment and reports correct coordinates
    Alignment loc = smithWaterman("TTTTGATTACATTTT", "CCGATTACACC");
    CHECK(loc.score == 14);
    CHECK(loc.alignedA == "GATTACA");
    CHECK(loc.aStart == 4 && loc.aEnd == 11 && loc.bStart == 2 && loc.bEnd == 9);

    // ungapped alignment strings are substrings at the reported coordinates
    std::string a = "ACACACTAGGGTTTACGT", b = "AGCACACAGGTTACGTCC";
    Alignment r = smithWaterman(a, b);
    CHECK(strip(r.alignedA) == a.substr(r.aStart, r.aEnd - r.aStart));
    CHECK(strip(r.alignedB) == b.substr(r.bStart, r.bEnd - r.bStart));
    CHECK(r.alignedA.size() == r.alignedB.size() && r.alignedA.size() == r.midline.size());

    if (failures == 0) std::cout << "All Smith-Waterman tests passed\n";
    return failures == 0 ? 0 : 1;
}
