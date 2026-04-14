#include <cassert>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

int main()
{
	unsigned long n = 12;
	int totalPerms = 0;
	for (unsigned long i=2; i<n; ++i) {
		for (unsigned long j=1; j<i; ++j) {
			std::vector<std::vector<int>> permutations;
				
			std::vector<int> indices(i);
			std::iota(indices.begin(), indices.end(), 0);
			
			std::vector<bool> select(i, false);
			std::fill(select.begin(), select.begin() + j, true);
			
			unsigned long fuck = 0;
			do {
				std::vector<int> chosen;
				for (unsigned long k = 0; k < i; ++k) {
					if (select[k]) {
						chosen.push_back(indices[k]);
					}
				}
				
				std::sort(chosen.begin(), chosen.end());
				
				do {
					permutations.push_back(chosen);
				} while (std::next_permutation(chosen.begin(), chosen.end()));

				fuck = fuck + 1;
			
			} while (std::prev_permutation(select.begin(), select.end()));
			
			totalPerms = permutations.size();
			std::cout << "n: " << i << " m: " << j << " Total Permutations: " << totalPerms << " Number of times code ran: " << fuck << std::endl;
			
			//for (int k = 0; k < totalPerms; ++k) {
			//	for (int l = 0; l < permutations[k].size(); ++l) {
			//		std::cout << permutations[k][l] << " ";
			//	}
			//	std::cout << std::endl;
			//}
		}
	}
}
