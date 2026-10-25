// smith_waterman.h -- local sequence alignment via dynamic programming + backtrace.
#pragma once
#include <cstddef>
#include <iostream>
#include <string>

/// Scoring parameters (linear gap penalty; gap should be negative).
struct Scoring {
    int match = 2;
    int mismatch = -1;
    int gap = -2;
};

/// Result of a local alignment. Start indices are 0-based inclusive, end indices exclusive.
struct Alignment {
    int score = 0;
    std::string alignedA, alignedB;  // gapped strings ('-' = gap)
    std::string midline;             // '|' for match, ' ' otherwise
    std::size_t aStart = 0, aEnd = 0, bStart = 0, bEnd = 0;
};

/// Smith-Waterman local alignment of [a] (rows) and [b] (columns).
/// H[i][j] = max(0, H[i-1][j-1]+s(a_i,b_j), H[i-1][j]+gap, H[i][j-1]+gap); the alignment is
/// backtraced from the highest-scoring cell until a cell with score 0 is reached.
/// Time and memory are O(|a|*|b|). Ties prefer diagonal, then gap in b, then gap in a.
Alignment smithWaterman(const std::string& a, const std::string& b, Scoring sc = Scoring{});

/// Pretty-prints an alignment in blocks of [width] columns, with 1-based coordinates.
void printAlignment(const Alignment& al, std::ostream& os = std::cout, std::size_t width = 60);
