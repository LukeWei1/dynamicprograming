// Runs Smith-Waterman: testAgainst vs targetGenome. Build with `make align`.
#include <algorithm>
#include <iostream>
#include "dna_data.h"
#include "smith_waterman.h"

int main() {
    Alignment al = smithWaterman(targetGenome, testAgainst);  // A = targetGenome, B = testAgainst
    printAlignment(al);
    auto pos = targetGenome.find(genomeSnippet);
    if (pos != std::string::npos)
        std::cout << "genomeSnippet occurs in targetGenome at " << pos + 1 << "-" << pos + genomeSnippet.size() << "\n";
    std::size_t matches = std::count(al.midline.begin(), al.midline.end(), '|');
    std::cout << "Identity: " << matches << "/" << al.midline.size() << "\n";
}
