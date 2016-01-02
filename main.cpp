//#include "SimpleObject.cpp"
#include <unistd.h>
#include <cstdio>
#include <sys/time.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>

using namespace std;

const int N = 9000; 			//number of triangles after the simplification
const double EPS = 1e-28;
const int LEFT = 47;
const int RIGHT = 12;

int sign(double x) {
	if(x < EPS && x > -EPS)
		return 0;
	if(x >= EPS)
		return 1;
	return -1;
}

struct Mat {
	double m[4][4];
	Mat() {
		for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
			m[i][j] = 0.;
		}
	}

	void plus(double a[4][4]) {
		for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
			m[i][j] += a[i][j];
		}
	}

	void clear() {
		  for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
			    m[i][j] = 0.;
		  }
	}

	double get(int i, int j) {
		return m[i][j];
	}

	void set(int i, int j, double x) {
		m[i][j] = x;
		return;
	}
};

class Point;

struct Tri {
	int id[3];
	Point* p[3];
	Tri* next;
	Tri* prev;
	Tri(int i0, int i1, int i2, Point* p0, Point* p1, Point* p2) {
		id[0] = i0;
		id[1] = i1;
		id[2] = i2;
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
	}
};

struct Pair {
	Point* p1;
	Point* p2;
	double cost;
	Pair* next;
	Pair* prev;
	Point* opt;
	Pair(Point* p1, Point* p2) {
		this->p1 = p1;
		this->p2 = p2;
		cost = 0.;
	}
	bool isEqual(Pair p) { 
		if(p1 == p.p1 && p2 == p.p2)
			return true;
		if(p1 == p.p2 && p2 == p.p1)
			return true;
		return false;
	}
};


bool cmp(Pair* pair1, Pair* pair2) {
	  return pair1->cost < pair2->cost;
}


class Point {
	public:
		double x,y,z;
		Mat mat;
		int index;
		Point* next;
		Point* prev;
		vector<Pair*> pairAtVertex;
		vector<Point*> neighbor;
		vector<Tri*> triAtVertex;
		Point() {}
		Point(double ix, double iy, double iz): x(ix), y(iy), z(iz)  {}

		Point operator+(const Point& p) const {
			return Point(this->x + p.x, this->y + p.y, this->z + p.z);
		}
		Point operator-(const Point& p) const {
			return Point(this->x - p.x, this->y - p.y, this->z - p.z);
		}
		Point operator*(const double c) const {
			return Point(this->x * c, this->y * c, this->z * c);
		}
		Point operator/(const double c) const {
			return Point(this->x / c, this->y / c, this->z / c);
		}

		double length() {
			return sqrt(x * x + y * y + z * z);
		}
		void normalize() {
			double len = length();
			x = x / len;
			y = y / len;
			z = z / len;
		}
		double distance(Point p) {
			return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y) + (z - p.z) * (z - p.z));
		}
};

class Model {
	public:
		int numberTri;
		Point* vertexHead;
		Tri* triangleHead;
		Pair* pairRoot;
		vector<Point*> v;
		vector<Pair*> pairs;			// store the valid pair
		vector<Tri*> triangle;			

		Model() {
			vertexHead = NULL;
			triangleHead = NULL;
			pairRoot = NULL;
			numberTri = 0;
		}

		void readObj() {
			//	ifstream fin("../../obj/bunny.fine.obj"); //	ifstream fin("../../obj/dragon.obj"); //	ifstream fin("cube.obj");
//			ifstream fin("../obj/block.obj");
//			ifstream fin("../obj/horse.fine.90k.obj");
			ifstream fin("../obj/Arma.obj");
			char c;
			double x, y, z;
			long int p[3];
			string temp;
			while (fin >> c) {
				if (c == 'v') {
					fin >> x >> y >> z;
					Point* pt = new Point(x,y,z);
					v.push_back(pt);
					pt->index = (int)v.size();
					pt->next = vertexHead;
					if(vertexHead != NULL)
						vertexHead->prev = pt;
					pt->prev = NULL;
					vertexHead = pt;
				} else if(c == 'f') {
					fin >> p[0] >> p[1] >> p[2];

					for(int i = 0; i < 3; i++)
						p[i]--;

					Tri* tri = new Tri(p[0], p[1], p[2], v[p[0]], v[p[1]], v[p[2]]);
					for(int i = 0; i < 3; i++) {
						tri->p[i]->triAtVertex.push_back(tri);
						for(int j = 0; j < 3; j++) {
							if(j == i)
								continue;
							bool isIn = false;
							for(int k = 0; k < (int)tri->p[i]->neighbor.size(); k++) {
								if(tri->p[i]->neighbor[k] == tri->p[j])
									isIn = true;
							}
							if(!isIn)
								tri->p[i]->neighbor.push_back(tri->p[j]);
						}
					}
					tri->next = triangleHead;
					if(triangleHead != NULL)
						triangleHead->prev = tri;
					tri->prev = NULL;
					triangleHead = tri;
					numberTri++;
				} else if(c == '#') {
					getline(fin, temp);
				}
			}
//			cout << " number of vertex is " << v.size() << endl;
			return;
		}

		void writeObj() {
//			ofstream fout("../../Light Tracing/src/DragonSimple/cubesimplified.obj");
			ofstream fout("../../Light Tracing/src/DragonSimple/Armasimplified.obj");
//			ofstream fout("../../Light Tracing/src/DragonSimple/blocksimplified3.obj");
//			ofstream fout("../../Light Tracing/src/DragonSimple/Horsesimplified.obj");
			Point* cur = vertexHead;
			int count = 1;
			while(cur != NULL) {
				fout << "v " << cur->x << ' ' << cur->y << ' ' << cur->z << endl;
				cur->index = count;
				cur = cur->next;
				count++;
			}
			Tri* tri = triangleHead;
			while(tri != NULL) {
				fout << "f " << tri->p[0]->index << ' ' << tri->p[1]->index << ' ' << tri->p[2]->index << endl;
				tri = tri->next;
			}
		}

		void init() { 
			Tri* tri = triangleHead;
			while(tri != NULL) {
				double ppT[4][4];
				for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) {
					ppT[i][j] = 0.;
				}
				computePPT(tri, ppT);
				for(int j = 0; j < 3; j++) {
					tri->p[j]->mat.plus(ppT);
				}
				tri = tri->next;
			}

			for(int i = 0; i < (int)v.size(); i++) {
				for(int j = 0; j < (int)(v[i]->neighbor.size()); j++) {
					if(v[i]->neighbor[j]->index > v[i]->index) {
						Pair* pr = new Pair(v[i], v[i]->neighbor[j]);
						pairs.push_back(pr);
						v[i]->pairAtVertex.push_back(pr);
						v[i]->neighbor[j]->pairAtVertex.push_back(pr);
					}
				}
			}
			return;
		}

		void computeCost() {
			for(size_t i = 0; i < pairs.size(); i++) {
				Point* p1 = pairs[i]->p1;
				Point* p2 = pairs[i]->p2;
				double QBar[4][4];
				for(int j = 0; j < 3; j++) for(int k = 0; k < 4; k++) {
					QBar[j][k] = p1->mat.get(j,k) + p2->mat.get(j,k);
				}
				for(int k = 0; k < 3; k++)
					QBar[3][k] = 0.;
				QBar[3][3] = 1.;
				double vBar[4] = {0.};
				double b[4] = {0,0,0,1};
				if(!solve(QBar, b, vBar)) {
					cout << "cannot solve!" << endl;
					vBar[0] = (p1->x + p2->x) / 2.;
					vBar[1] = (p1->y + p2->y) / 2.;
					vBar[2] = (p1->z + p2->z) / 2.;
					vBar[3] = 1.;
				}
				pairs[i]->opt = new Point(vBar[0], vBar[1], vBar[2]);
				for(int j = 0; j < 4; j++) for(int k = 0; k < 4; k++) {
					QBar[j][k] = p1->mat.get(j,k) + p2->mat.get(j,k);
				}
				pairs[i]->cost = computeQuadratic(QBar, vBar);
			}
		}

		void buildList() {
			pairRoot = pairs[0];
			pairs[0]->prev = NULL;
			for(size_t i = 0; i < pairs.size() - 1; i++) {
				pairs[i]->next = pairs[i + 1];
				pairs[i + 1]->prev = pairs[i];
			}
			pairs[pairs.size() - 1]->next = NULL;
		}

		void simplify(int targetNumber) {
			int numberOfVertex = (int)v.size();
			while(numberOfVertex > targetNumber) {
				  if(numberOfVertex < 6000) {
					    simplifyOne(targetNumber, numberOfVertex);
					    break;
				  }
				  vector<Pair*> sortPairs;
				  Pair* cur = pairRoot;
				  while(cur != NULL) {
					    sortPairs.push_back(cur);
					    cur = cur->next;
				  }
				  sort(sortPairs.begin(),sortPairs.end(),cmp);
				  Pair* minCostPair = sortPairs[0];
				  double minCost = minCostPair->cost;
				  double maxCost;
				  if(numberOfVertex > (int)(v.size() * 0.6)) {
					    maxCost = minCost * 18.0;
				  } else if(numberOfVertex > (int)(v.size() * 0.4)) {
					    maxCost = minCost * 6.2;
				  } else if(numberOfVertex > (int)(v.size() * 0.2)) {
					    maxCost = minCost * 4.4;
				  } else {
					    maxCost = minCost * 2.5;
				  }
//				  maxCost = minCost * 2.2;
				  vector<bool> conflictVertex;
				  for(int i = 0; i < (int)v.size(); i++) {
					    conflictVertex.push_back(false);
				  }
				  int index = 0;
				  while(minCostPair->cost < maxCost && index < (int)sortPairs.size() * 0.8) {
					    if(minCostPair->cost > maxCost)
							break;
//					    cout << numberOfVertex << endl;
					    if(conflictVertex[minCostPair->p1->index - 1] || conflictVertex[minCostPair->p2->index - 1]) {
							index++;
							minCostPair = sortPairs[index];
							continue;
					    }

					    conflictVertex[minCostPair->p1->index - 1] = true;
					    conflictVertex[minCostPair->p2->index - 1] = true;

					    Point* p1 = minCostPair->p1;
					    Point* p2 = minCostPair->p2; //we will delete p2 and leave p1
					    Point* pv = minCostPair->opt;
					    Tri* tri1 = NULL;


					    //delete the pair p1p2 from the pair list of p1
					    for (int i = 0; i < (int)p1->pairAtVertex.size(); i++) {
							if(p1->pairAtVertex[i]->p1 == p2 || p1->pairAtVertex[i]->p2 == p2) {
								  p1->pairAtVertex.erase(p1->pairAtVertex.begin() + i);
								  i--;
							}
					    }


					    //changing the position of p1
					    p1->x = pv->x;
					    p1->y = pv->y;
					    p1->z = pv->z;

					    //for each p3 update 
					    for(int j = 0; j < (int)p1->triAtVertex.size(); j++) {
							bool del = false;
							for(int k = 0; k < 3; k++) {
								  if(p1->triAtVertex[j]->p[k] == p2) {

									    tri1 = p1->triAtVertex[j];

									    //delete tri1 from tri list of p1
									    for (int i = 0; i < (int)p1->triAtVertex.size(); i++) {
											if(p1->triAtVertex[i] == tri1) {
												  p1->triAtVertex.erase(p1->triAtVertex.begin() + i);
												  del = true;
											}
									    }

									    //find p3
									    Point* p3;
									    for (int i = 0; i < 3; i++) {
											if(tri1->p[i] != p1 && tri1->p[i] != p2)
												  p3 = tri1->p[i];
									    }

									    //modify the pair list of p3 
									    for (int i = 0; i < (int)p3->pairAtVertex.size(); i++) {
											if(p3->pairAtVertex[i]->p1 == p2 || p3->pairAtVertex[i]->p2 == p2) {
												  p3->pairAtVertex.erase(p3->pairAtVertex.begin() + i);
												  i--;
											}
									    }


									    //modify the tri list of p3 
									    for (int i = 0; i < (int)p3->triAtVertex.size(); i++) {
											if(p3->triAtVertex[i] == tri1) {
												  p3->triAtVertex.erase(p3->triAtVertex.begin() + i);
												  i--;
											}
									    }

									    //remove the tri
									    if(tri1 == triangleHead)
											triangleHead = tri1->next;
									    if(tri1->next != NULL) {
											tri1->next->prev = tri1->prev;
									    }
									    if(tri1->prev != NULL) {
											tri1->prev->next = tri1->next;
									    }

									    //update the pair list of p2, delete p2p3 from it and delete this pair
									    for(int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
											Pair* dp;
											if(p2->pairAtVertex[i]->p1 == p3 || p2->pairAtVertex[i]->p2 == p3) {
												  dp = p2->pairAtVertex[i];
												  if(dp == pairRoot)
													    pairRoot = dp->next;
												  if(dp->next != NULL) {
													    dp->next->prev = dp->prev;
												  }
												  if(dp->prev != NULL) {
													    dp->prev->next = dp->next;
												  }
												  p2->pairAtVertex.erase(p2->pairAtVertex.begin() + i);
												  i--;
												  //									cout << "delete pair " << dp->p1->index << ' ' << dp->p2->index << endl;
											}
									    }
								  }
							}
							if(del)
								  j--;
					    }

					    //delete pair p1p2
					    for(int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
							Pair* dp;
							if(p2->pairAtVertex[i]->p1 == p1 || p2->pairAtVertex[i]->p2 == p1) {
								  dp = p2->pairAtVertex[i];
								  if(dp == pairRoot)
									    pairRoot = dp->next;
								  if(dp->next != NULL) {
									    dp->next->prev = dp->prev;
								  }
								  if(dp->prev != NULL) {
									    dp->prev->next = dp->next;
								  }
								  p2->pairAtVertex.erase(p2->pairAtVertex.begin() + i);
								  i--;
								  //						cout << "delete primal pair " << dp->p1->index << ' ' << dp->p2->index << endl;
							}
					    }

					    //modify the pair list of p1 and the pairs that will be attached to p1
					    for (int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
							Point* p2Neighbor;
							if(p2->pairAtVertex[i]->p1 == p2) {
								  p2->pairAtVertex[i]->p1 = p1;
								  p2Neighbor = p2->pairAtVertex[i]->p2;
							} else {
								  p2->pairAtVertex[i]->p2 = p1;
								  p2Neighbor = p2->pairAtVertex[i]->p1;
							}
							p1->pairAtVertex.push_back(p2->pairAtVertex[i]);
							for (int j = 0; j < (int)p2Neighbor->pairAtVertex.size(); j++) {
								  if(p2Neighbor->pairAtVertex[j]->p1 == p2) {
									    p2Neighbor->pairAtVertex[j]->p1 = p1;
								  }
								  if(p2Neighbor->pairAtVertex[j]->p2 == p2) {
									    p2Neighbor->pairAtVertex[j]->p2 = p1;
								  }
							}

							for(size_t j = 0; j < p2Neighbor->pairAtVertex.size(); j++) {
								  for(size_t k = j + 1; k < p2Neighbor->pairAtVertex.size(); k++) {
									    if(p2Neighbor->pairAtVertex[j]->p1 == p2Neighbor->pairAtVertex[k]->p1 && p2Neighbor->pairAtVertex[j]->p2 == p2Neighbor->pairAtVertex[k]->p2) {
											p2Neighbor->pairAtVertex.erase(p2Neighbor->pairAtVertex.begin() + k);
											k--;
									    }
								  }
							}

					    }


					    for(int i = 0; i < (int)p2->triAtVertex.size(); i++) {
							bool isP3Tri = false;
							for(int j = 0; j < 3; j++) {
								  if(p2->triAtVertex[i]->p[j] == p1)
									    isP3Tri = true;
							}
							if(isP3Tri)
								  continue;
							if(p2->triAtVertex[i]->p[0] == p2)
								  p2->triAtVertex[i]->p[0] = p1;
							if(p2->triAtVertex[i]->p[1] == p2)
								  p2->triAtVertex[i]->p[1] = p1;
							if(p2->triAtVertex[i]->p[2] == p2)
								  p2->triAtVertex[i]->p[2] = p1;
							p1->triAtVertex.push_back(p2->triAtVertex[i]);
					    }

					    if(p2 == vertexHead)
							vertexHead = p2->next;
					    if(p2->next != NULL) {
							p2->next->prev = p2->prev;
					    }
					    if(p2->prev != NULL) {
							p2->prev->next = p2->next;
					    }
					    numberOfVertex--;

					    p1->mat.clear();
					    for(int i = 0; i < (int)p1->triAtVertex.size(); i++) {
							double ppT[4][4];
							computePPT(p1->triAtVertex[i], ppT);
							p1->mat.plus(ppT);
					    }

					    for(size_t i = 0; i < p1->pairAtVertex.size(); i++) {
							Point* tp1 = p1->pairAtVertex[i]->p1;
							Point* tp2 = p1->pairAtVertex[i]->p2;
							double QBar[4][4];
							for(int j = 0; j < 3; j++) for(int k = 0; k < 4; k++) {
								  QBar[j][k] = tp1->mat.get(j,k) + tp2->mat.get(j,k);
							}
							for(int k = 0; k < 3; k++)
								  QBar[3][k] = 0.;
							QBar[3][3] = 1.;
							double vBar[4] = {0.};
							double b[4] = {0,0,0,1};
							if(!solve(QBar, b, vBar)) {
								  cout << "cannot solve!" << endl;
								  vBar[0] = (tp1->x + tp2->x) / 2.;
								  vBar[1] = (tp1->y + tp2->y) / 2.;
								  vBar[2] = (tp1->z + tp2->z) / 2.;
								  vBar[3] = 1.;
							}
							p1->pairAtVertex[i]->opt = new Point(vBar[0], vBar[1], vBar[2]);
							for(int j = 0; j < 4; j++) for(int k = 0; k < 4; k++) {
								  QBar[j][k] = tp1->mat.get(j,k) + tp2->mat.get(j,k);
							}
							p1->pairAtVertex[i]->cost = computeQuadratic(QBar, vBar);
					    }

					    for(size_t i = 0; i < p1->pairAtVertex.size(); i++) {
							for(size_t j = i + 1; j < p1->pairAtVertex.size(); j++) {
								  if(p1->pairAtVertex[i]->p1 == p1->pairAtVertex[j]->p1 && p1->pairAtVertex[i]->p2 == p1->pairAtVertex[j]->p2) {
									    p1->pairAtVertex.erase(p1->pairAtVertex.begin() + j);
									    j--;
								  }
							}
					    }

					    index++;
					    minCostPair = sortPairs[index];
				  }
				  cout << "vertex deleted " << index << endl;
			}
		}

		void simplifyOne(int NumberOfTri, int remain) {
			for (int round = 0; round < remain - NumberOfTri + 1; round++) {
				Pair* minCostPair = pairRoot;
				Pair* cur = pairRoot->next;
				while(cur != NULL) {
					if(sign(cur->cost - minCostPair->cost) < 0) {
						minCostPair = cur;
					}
					cur = cur->next;
				}

				Point* p1 = minCostPair->p1;
				Point* p2 = minCostPair->p2; //we will delete p2 and leave p1
				Point* pv = minCostPair->opt;
				Tri* tri1 = NULL;

				cout << "round " << round << endl;

				//delete the pair p1p2 from the pair list of p1
				for (int i = 0; i < (int)p1->pairAtVertex.size(); i++) {
					if(p1->pairAtVertex[i]->p1 == p2 || p1->pairAtVertex[i]->p2 == p2) {
						p1->pairAtVertex.erase(p1->pairAtVertex.begin() + i);
						i--;
					}
				}


				//changing the position of p1
				p1->x = pv->x;
				p1->y = pv->y;
				p1->z = pv->z;

				//for each p3 update 
				for(int j = 0; j < (int)p1->triAtVertex.size(); j++) {
					bool del = false;
					for(int k = 0; k < 3; k++) {
						if(p1->triAtVertex[j]->p[k] == p2) {

							tri1 = p1->triAtVertex[j];

							//delete tri1 from tri list of p1
							for (int i = 0; i < (int)p1->triAtVertex.size(); i++) {
								if(p1->triAtVertex[i] == tri1) {
//									cout << p1->triAtVertex[i]->p[0]->index << ' ' << p1->triAtVertex[i]->p[1]->index << ' ' << p1->triAtVertex[i]->p[2]->index << endl;
									p1->triAtVertex.erase(p1->triAtVertex.begin() + i);
									del = true;
								}
							}

							//find p3
							Point* p3;
							for (int i = 0; i < 3; i++) {
								if(tri1->p[i] != p1 && tri1->p[i] != p2)
									p3 = tri1->p[i];
							}
							//modify the pair list of p3 
							for (int i = 0; i < (int)p3->pairAtVertex.size(); i++) {
								if(p3->pairAtVertex[i]->p1 == p2 || p3->pairAtVertex[i]->p2 == p2) {
									p3->pairAtVertex.erase(p3->pairAtVertex.begin() + i);
									i--;
								}
							}

							//modify the tri list of p3 
							for (int i = 0; i < (int)p3->triAtVertex.size(); i++) {
								if(p3->triAtVertex[i] == tri1){
									p3->triAtVertex.erase(p3->triAtVertex.begin() + i);
									i--;
								}
							}

							//remove the tri
							if(tri1 == triangleHead)
								triangleHead = tri1->next;
							if(tri1->next != NULL) {
								tri1->next->prev = tri1->prev;
							}
							if(tri1->prev != NULL) {
								tri1->prev->next = tri1->next;
							}

							//update the pair list of p2, delete p2p3 from it and delete this pair
							for(int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
								Pair* dp;
								if(p2->pairAtVertex[i]->p1 == p3 || p2->pairAtVertex[i]->p2 == p3) {
									dp = p2->pairAtVertex[i];
									if(dp == pairRoot)
										pairRoot = dp->next;
									if(dp->next != NULL) {
										dp->next->prev = dp->prev;
									}
									if(dp->prev != NULL) {
										dp->prev->next = dp->next;
									}
									p2->pairAtVertex.erase(p2->pairAtVertex.begin() + i);
									i--;
//									cout << "delete pair " << dp->p1->index << ' ' << dp->p2->index << endl;
								}
							}
						}
					}
					if(del)
						  j--;
				}

				//delete pair p1p2
				for(int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
					Pair* dp;
					if(p2->pairAtVertex[i]->p1 == p1 || p2->pairAtVertex[i]->p2 == p1) {
						dp = p2->pairAtVertex[i];
						if(dp == pairRoot)
							pairRoot = dp->next;
						if(dp->next != NULL) {
							dp->next->prev = dp->prev;
						}
						if(dp->prev != NULL) {
							dp->prev->next = dp->next;
						}
						p2->pairAtVertex.erase(p2->pairAtVertex.begin() + i);
//						cout << "delete primal pair " << dp->p1->index << ' ' << dp->p2->index << endl;
						i--;
					}
				}

				//modify the pair list of p1 and the pairs that will be attached to p1
				for (int i = 0; i < (int)p2->pairAtVertex.size(); i++) {
					Point* p2Neighbor;
					if(p2->pairAtVertex[i]->p1 == p2) {
						p2->pairAtVertex[i]->p1 = p1;
						p2Neighbor = p2->pairAtVertex[i]->p2;
					} else {
						p2->pairAtVertex[i]->p2 = p1;
						p2Neighbor = p2->pairAtVertex[i]->p1;
					}
					p1->pairAtVertex.push_back(p2->pairAtVertex[i]);
					for (int j = 0; j < (int)p2Neighbor->pairAtVertex.size(); j++) {
						if(p2Neighbor->pairAtVertex[j]->p1 == p2) {
							p2Neighbor->pairAtVertex[j]->p1 = p1;
						}
						if(p2Neighbor->pairAtVertex[j]->p2 == p2) {
							p2Neighbor->pairAtVertex[j]->p2 = p1;
						}
					}

					for(size_t j = 0; j < p2Neighbor->pairAtVertex.size(); j++) {
						for(size_t k = j + 1; k < p2Neighbor->pairAtVertex.size(); k++) {
							if(p2Neighbor->pairAtVertex[j]->p1 == p2Neighbor->pairAtVertex[k]->p1 && p2Neighbor->pairAtVertex[j]->p2 == p2Neighbor->pairAtVertex[k]->p2) {
								p2Neighbor->pairAtVertex.erase(p2Neighbor->pairAtVertex.begin() + k);
								k--;
							}
						}
					}

				}


				for(int i = 0; i < (int)p2->triAtVertex.size(); i++) {
					bool isP3Tri = false;
					for(int j = 0; j < 3; j++) {
						if(p2->triAtVertex[i]->p[j] == p1)
							isP3Tri = true;
					}
					if(isP3Tri)
						continue;
					if(p2->triAtVertex[i]->p[0] == p2)
						p2->triAtVertex[i]->p[0] = p1;
					if(p2->triAtVertex[i]->p[1] == p2)
						p2->triAtVertex[i]->p[1] = p1;
					if(p2->triAtVertex[i]->p[2] == p2)
						p2->triAtVertex[i]->p[2] = p1;
					p1->triAtVertex.push_back(p2->triAtVertex[i]);
				}

				if(p2 == vertexHead)
					vertexHead = p2->next;
				if(p2->next != NULL) {
					p2->next->prev = p2->prev;
				}
				if(p2->prev != NULL) {
					p2->prev->next = p2->next;
				}

				p1->mat.clear();
				for(int i = 0; i < (int)p1->triAtVertex.size(); i++) {
					  double ppT[4][4];
					  computePPT(p1->triAtVertex[i], ppT);
					  p1->mat.plus(ppT);
				}

				for(size_t i = 0; i < p1->pairAtVertex.size(); i++) {
					  Point* tp1 = p1->pairAtVertex[i]->p1;
					  Point* tp2 = p1->pairAtVertex[i]->p2;
					  double QBar[4][4];
					  for(int j = 0; j < 3; j++) for(int k = 0; k < 4; k++) {
						    QBar[j][k] = tp1->mat.get(j,k) + tp2->mat.get(j,k);
					  }
					  for(int k = 0; k < 3; k++)
						    QBar[3][k] = 0.;
					  QBar[3][3] = 1.;
					  double vBar[4] = {0.};
					  double b[4] = {0,0,0,1};
					  if(!solve(QBar, b, vBar)) {
						    cout << "cannot solve!" << endl;
						    vBar[0] = (tp1->x + tp2->x) / 2.;
						    vBar[1] = (tp1->y + tp2->y) / 2.;
						    vBar[2] = (tp1->z + tp2->z) / 2.;
						    vBar[3] = 1.;
					  }
					  p1->pairAtVertex[i]->opt = new Point(vBar[0], vBar[1], vBar[2]);
					  for(int j = 0; j < 4; j++) for(int k = 0; k < 4; k++) {
						    QBar[j][k] = tp1->mat.get(j,k) + tp2->mat.get(j,k);
					  }
					  p1->pairAtVertex[i]->cost = computeQuadratic(QBar, vBar);
				}

				for(size_t i = 0; i < p1->pairAtVertex.size(); i++) {
					for(size_t j = i + 1; j < p1->pairAtVertex.size(); j++) {
						if(p1->pairAtVertex[i]->p1 == p1->pairAtVertex[j]->p1 && p1->pairAtVertex[i]->p2 == p1->pairAtVertex[j]->p2) {
							p1->pairAtVertex.erase(p1->pairAtVertex.begin() + j);
							j--;
						}
					}
				}
			}
		}



		void computePPT(Tri* tri, double a[4][4]) { 
			Point norm = prod(*v[tri->id[1]] - *v[tri->id[0]], *v[tri->id[2]] - *v[tri->id[0]]);
			norm.normalize();
			double p[4] = {norm.x, norm.y, norm.z, dot(norm, *(tri->p[0])) * (-1)};
			for(int i = 0; i < 4; i++) for(int j = 0; j < 4; j++) { 
				a[i][j] = p[i] * p[j];
			}
		}

		bool solve(double A[4][4], double B[4], double x[4]) {
			double a[4][4];

			double b[4];
			for(int i = 0; i < 4; i++) {
				b[i] = B[i];
				for(int j = 0; j < 4; j++) { 
					a[i][j] = A[i][j];
				}
			}
			for(int k = 0; k < 4; k++) {
				if(sign(a[k][k]) == 0)
					return false;
				for(int i = 0; i < 4; i++) {
					if(i == k)
						continue;
					double c = a[i][k] / a[k][k];
					for(int j = k + 1; j < 4; j++) { 
						a[i][j] = a[i][j] - c * a[k][j];
					}
					b[i] = b[i] - c * b[k];
				}
			}
			for(int i = 0; i < 4; i++) {
				x[i] = b[i] / a[i][i];
			}
			return true;
		}

		double computeQuadratic(double A[4][4], double x[4]) {
			double res = 0.;
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < 4; j++) {
					res += x[i] * A[i][j] * x[j];
				}
			}
			return res;
		}

		Point prod(Point p1, Point p2) {
			return Point(p1.y * p2.z - p2.y * p1.z, p1.z * p2.x - p2.z * p1.x, p1.x * p2.y - p2.x * p1.y);
		}

		double dot(Point p1, Point p2) {
			return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
		}

};

int main() {
	  struct timeval tpstart,tpend;
	  double timeuse;

	  gettimeofday(&tpstart,NULL);

	  Model m = Model();
	  m.readObj();
	  m.init();
	  m.computeCost();
	  //	m.buildHeap();
	  m.buildList();
	  m.simplify(N);
	  m.writeObj();
	  cout << "<======= Simplification finished =======>" << endl;

	  gettimeofday(&tpend,NULL);

	  timeuse=1000000*(tpend.tv_sec-tpstart.tv_sec)+tpend.tv_usec-tpstart.tv_usec;

	  timeuse/=1000000;

	  printf("time used is %lf s\n",timeuse);
	  return 0;
}
