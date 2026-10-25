// smith_waterman.cpp -- implementation (see smith_waterman.h for documentation).
#include "smith_waterman.h"

#include <algorithm>
#include <vector>

// Takes the two sequences and the scoring settings, returns an Alignment.
Alignment smithWaterman(const std::string& a, const std::string& b, Scoring sc) {
    // m, n = sequence lengths; W = table row width (n + 1 for the empty-prefix column).
    const std::size_t m = a.size(), n = b.size(), W = n + 1;

    // Flat (m+1) x (n+1) score table of zeros; cell (i, j) is H[i * W + j].
    std::vector<int> H((m + 1) * W, 0);

    // Highest score so far and where it is.
    int best = 0;
    std::size_t bi = 0, bj = 0;

    // Fill the table. Table index i matches a[i-1], and j matches b[j-1].
    for (std::size_t i = 1; i <= m; ++i) {
        for (std::size_t j = 1; j <= n; ++j) {
            // Align a[i-1] with b[j-1]: match or mismatch score.
            int diag = H[(i - 1) * W + j - 1] + (a[i - 1] == b[j - 1] ? sc.match : sc.mismatch);
            // Gap in b (move down a row).
            int up = H[(i - 1) * W + j] + sc.gap;
            // Gap in a (move right a column).
            int left = H[i * W + j - 1] + sc.gap;

            // The 0 restarts the alignment, which is what makes it local.
            int v = std::max({0, diag, up, left});
            H[i * W + j] = v;

            // Keep the first maximum found.
            if (v > best) { best = v; bi = i; bj = j; }
        }
    }

    // The alignment ends at the highest-scoring cell.
    Alignment res;
    res.score = best;
    res.aEnd = bi;
    res.bEnd = bj;

    // Traceback from the max cell until the score hits 0.
    std::size_t i = bi, j = bj;
    while (i > 0 && j > 0 && H[i * W + j] > 0) {
        int cur = H[i * W + j];
        int s = (a[i - 1] == b[j - 1]) ? sc.match : sc.mismatch;

        // Work out which neighbour gave this score (ties prefer the diagonal).
        if (cur == H[(i - 1) * W + j - 1] + s) {
            // Diagonal: aligned pair; '|' marks a match.
            res.alignedA += a[i - 1];
            res.alignedB += b[j - 1];
            res.midline += (a[i - 1] == b[j - 1]) ? '|' : ' ';
            --i; --j;
        } else if (cur == H[(i - 1) * W + j] + sc.gap) {
            // Up: gap in B.
            res.alignedA += a[i - 1];
            res.alignedB += '-';
            res.midline += ' ';
            --i;
        } else {
            // Left: gap in A.
            res.alignedA += '-';
            res.alignedB += b[j - 1];
            res.midline += ' ';
            --j;
        }
    }

    // Strings were built backward, so reverse them.
    std::reverse(res.alignedA.begin(), res.alignedA.end());
    std::reverse(res.alignedB.begin(), res.alignedB.end());
    std::reverse(res.midline.begin(), res.midline.end());

    // (i, j) is now the 0-based start of the alignment.
    res.aStart = i;
    res.bStart = j;
    return res;
}
//Prints the header. Start positions are converted from 0-based to 1-based. The end is already exclusive, 
//which is the same number as a 1-based inclusive end.
void printAlignment(const Alignment& al, std::ostream& os, std::size_t width) {
    os << "Score: " << al.score << "\n"
       << "Sequence A: " << al.aStart + 1 << "-" << al.aEnd << "\n"
       << "Sequence B: " << al.bStart + 1 << "-" << al.bEnd << "\n\n";
    std::size_t ia = al.aStart + 1, ib = al.bStart + 1;
    for (std::size_t p = 0; p < al.alignedA.size(); p += width) {
        std::string sa = al.alignedA.substr(p, width), sb = al.alignedB.substr(p, width);
        std::size_t na = sa.size() - std::count(sa.begin(), sa.end(), '-');
        std::size_t nb = sb.size() - std::count(sb.begin(), sb.end(), '-');
        os << "A " << ia << "\t" << sa << "\n  \t" << al.midline.substr(p, width) << "\nB " << ib << "\t" << sb << "\n\n";
        ia += na; ib += nb;
    }
}
