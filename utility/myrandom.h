#include <iterator>
#include <set>

float fran(float mult = 1.0);
double dran(double mult = 1.0);
int iran(int max);
unsigned int plran(float l, float u, float r);

template <class T>
T getRandomElement(std::set<T> S) {
	unsigned int pidx = iran(S.size());
	typename std::set<T>::iterator iter = S.begin();
	for (unsigned int k(0); k < pidx; ++k) {
		++iter;
	}
	return(*iter);
}
