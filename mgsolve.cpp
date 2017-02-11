#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#define _USE_MATH_DEFINES

using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::ofstream;

//class grid: used to manange the operation of matrix
class grid {
private:
	const size_t l;
	const long height;
	const long width;
	vector<double> u;

	vector<double> get_u() const
	{	return u;
	}

public:
	//Constructor with level
	grid (const size_t& i = 1) 
		: l (i) , height (pow (2 , i) + 1) , width (pow (2 , i) + 1)
	{	vector<double> temp(height * width , 0);
		u = temp;
	}

	//Constructor with width and height
	grid (const long& i , const long& j) : l (0) , height (i) , width (j)
	{	vector<double> temp(height * width , 0);
		u = temp;
	}

	//Copy constructor
	grid (const grid &source) : l (source.get_level()) , height (source.get_height()) , width (source.get_width())
	{	u = source.get_u();
	}
	
	~grid()
	{}

	//Return value at node [row , col]
	double get (long row , long col)
	{	double res = 0;
		if (row >= height || col >= width)
		{	cout << "Error: out of boundary!" << endl;
			return -1;
		}
		else
		{	vector<double>::iterator ite_u = u.begin();
			ite_u += row * width + col;
			res = *ite_u;
			return res;
		}
	}

	//Edit node [row , col] with val
	void set (long row , long col , double val)
	{	if (row >= height || col >= width)
		{	cout << "Error: out of boundary!" << endl;
		}
		else
		{	vector<double>::iterator ite_u = u.begin();
			ite_u += row * width + col;
			*ite_u = val;
		}
	}

	//Return the level of the grid
	size_t get_level () const
	{	return l;
	}

	//return the number of nodes of the grid
	long get_node() const
	{	if (width == pow (2 , l) + 1 && height == width)
		{	return width;
		} else
		{	cout << "Inappropriate call of get_node()!" << endl;
			return 0;
		}
	}

	//return the width of grid
	long get_width() const
	{	return width;
	}

	//return the height of grid
	long get_height() const
	{	return height;
	}

	//Print the gird on the screen and output it to file "output.txt"
	void output()
	{	/*for (long i = 0 ; i < width; ++i)
		{	for (long j = 0 ; j < height; ++j)
			{	cout << get (i , j) << " ";
			}
			cout << endl;
		}*/

		ofstream out("solution.txt");

		if(!out.is_open()){
			cerr << "Error: Failed to open solution.txt!" << endl;
			return ;
		}

		double h = 1 / pow (2 , l);
		for(long i = 0 ; i < width ; ++i){
			for(long j = 0 ; j < height ; ++j)
			{	out << j * h << " " << i * h << " " << get(i , j) << endl;
			}
		}		

		out.close();
	}

	//set all nodes to 0
	void zero()
	{	for (double i = 0 ; i < width; ++i)
		{	for (double j = 0 ; j < height; ++j)
			{	set(i,j,0);
			}
		}
	}

	//copy assignment
	grid& operator = (const grid& source)
	{	if (l == source.get_level() && width == source.get_width() && height == source.get_height())
		{	u = source.get_u();
		}
		
		return *this;
	}

	//Return the quadratic sum of all nodes
	double quad_sum()
	{	double res = 0;

		for (long i = 0 ; i < width ; ++i)
		{	for (long j = 0 ; j < height ; ++j)
			{	res += get(i , j) * get(i , j);
			}
		}
		return res;
	}
};

//Return the difference: a-b
grid diff(grid a , grid b)
{	if (a.get_width() != b.get_width() || a.get_height()!= b.get_height())
	{	cout << "Error: Require grids of same size!" << endl;
		return a;
	}
	
	grid e(a);
	for (long i = 0 ; i < a.get_width() ; ++i)
	{	for (long j = 0 ; j < a.get_height() ; ++j)
		e.set(i , j , a.get(i,j) - b.get(i,j));
	}

	return e;
};

//class mgsolve : multigrid solver
class mgsolve {
private:
	const size_t l;
	size_t n;
	vector <grid> u;
	vector <grid> r;
	vector <grid> e;
	size_t current_level;

	//Restriction from rh to r2h
	void restriction(vector <grid>::iterator rh , vector <grid>::iterator r2h)
	{	double temp = 0;

		r2h -> get_level();
		for (long i = 1 ; i < r2h->get_node() - 1 ; ++i)
		{	for (long j = 1 ; j < r2h->get_node() - 1 ; ++j)
			{	temp = ( rh->get((2 * i + 1), (2 * j + 1)) +
						 rh->get((2 * i + 1), (2 * j    )) * 2 +
						 rh->get((2 * i + 1), (2 * j - 1)) +
						 rh->get((2 * i    ), (2 * j + 1)) * 2 +
						 rh->get((2 * i    ), (2 * j    )) * 4 +
						 rh->get((2 * i    ), (2 * j - 1)) * 2 +
						 rh->get((2 * i - 1), (2 * j + 1)) +
						 rh->get((2 * i - 1), (2 * j    )) * 2 +
						 rh->get((2 * i - 1), (2 * j - 1))
						) / 16;
				r2h->set (i , j , temp);
			}
		}
	}

	//Interpolation from r2h to rh
	void interpolation (vector <grid>::iterator r2h , vector <grid>::iterator rh)
	{	double temp = 0;

		for (int i = 0 ; i < r2h->get_node() - 1 ; ++i)
		{	for (int j = 0 ; j < r2h->get_node() - 1 ; ++j)
			{	temp = r2h->get (i , j);
				rh->set (2 * i , 2 * j , temp);
				temp = (r2h->get (i , j) + r2h->get (i + 1 , j)) / 2;
				rh->set (2 * i + 1 , 2 * j , temp);
				temp = (r2h->get (i , j) + r2h->get (i , j + 1)) / 2;
				rh->set (2 * i , 2 * j + 1 , temp);
				temp = (r2h->get (i , j) + r2h->get (i + 1 , j) + r2h->get (i , j + 1) + r2h->get (i + 1 , j + 1)) / 4;
				rh->set (2 * i + 1 , 2 * j + 1 , temp);
			}
		}
	}

	//Red-black Gauss-Seidel Algorithm for -delat uh = fh
	void gs (vector <grid>::iterator uh , vector <grid>::iterator fh)
	{	long node = uh->get_node() - 1;
		double temp = 0;

		for (int i = 1 ; i < node ; ++i)
		{	for (int j = 1 ; j < node ; ++j)
			{	if ((i + j) % 2 == 0)	// red point
				{	temp =	(fh->get(i , j) / (node * node) + 
							 uh->get (i - 1 , j) + 
							 uh->get (i + 1 , j) + 
							 uh->get (i , j - 1) + 
							 uh->get (i , j + 1)) / 4;
					uh->set(i , j , temp);
				}
			}
		}	

		for (int i = 1 ; i < node ; ++i)
		{	for (int j = 1 ; j < node ; ++j)
			{	if ((i + j) % 2 == 1)	// black point
				{	temp =	(fh->get(i , j) / (node * node) + 
							 uh->get (i - 1 , j) + 
							 uh->get (i + 1 , j) + 
							 uh->get (i , j - 1) + 
							 uh->get (i , j + 1)) / 4;
					uh->set(i , j , temp);
				}
			}
		}
	}

	//Compute the residual: r = f - (-delat u)
	void compute_residual (vector<grid>::iterator fh , vector<grid>::iterator uh , vector<grid>::iterator rh)
	{	double temp = 0;
		long node = rh->get_node() - 1;

		for (int i = 1 ; i < node ; ++i)
		{	for (int j = 1 ; j < node ; ++j)
			{	temp = fh->get(i , j) + 
						(uh->get(i , j - 1) + uh->get(i , j + 1) + 
						 uh->get(i - 1 , j) + uh->get(i + 1 , j) -
						 uh->get(i , j) * 4) * 
						 node * node;
				rh->set (i , j , temp);
			}
		}
	}

	//Correct the approximation uh with the error eh
	void correct (vector<grid>::iterator uh , vector<grid>::iterator eh)
	{	double temp = 0;

		for (long i = 0 ; i < uh->get_node() ; ++i )
		{	for (long j = 0 ; j < uh->get_node() ; ++j )
			{	temp = uh->get(i , j) + eh->get(i , j);
				uh->set(i , j , temp);
				eh->set(i , j , 0);
			}
		}
	}

	//Multigrid algorithm
	void mgm (vector<grid>::iterator uh , vector<grid>::iterator fh, size_t v1 = 2 , size_t v2 = 1)
	{	for (size_t i = 0 ; i < v1 ; ++i)
		{	gs(uh , fh);
		}

		vector<grid>::iterator rh = r.begin() + (l - uh -> get_level());
		compute_residual(fh , uh , rh);
		
		vector<grid>::iterator f2h = fh + 1;
		restriction(rh , f2h);
		restriction(uh , uh+1);

		vector<grid>::iterator eh = e.begin() + (l - uh -> get_level());
		vector<grid>::iterator e2h = eh + 1;
		
		if (uh-> get_level() == 2)
		{	gs (e2h , f2h);
		} else
		{	e2h->zero();
			mgm(e2h , f2h , 2 , 1);
		}

		interpolation (e2h , eh);
		correct (uh , eh);
		
		for (size_t i = 0 ; i < v2 ; ++i)
		{	gs(uh , fh);
		}
	}

public:
	//Constructor with level and number of cycles
	mgsolve (const size_t& x = 2 , const size_t& y = 1) : l (x) , n (y)
	{	u.push_back(*new grid(l));
		r.push_back(*new grid(l));
		e.push_back(*new grid(l));
		
		current_level = 0;
		double h = 1 / pow (2 , l);
		
		//Dirichlet BCs
		for (long i = 0 ; i < pow (2 , l) + 1 ; ++i)
		{	for (long j = 0 ; j < pow (2 , l) + 1 ; ++j)
			{	if ( i == 0 || i == pow (2 , l) || j == 0 || j == pow (2 , l))
				{	double temp = sin(M_PI * h * i) * sinh (M_PI * h * j);
					u.begin()->set (i , j , temp);
				}
			}
		}

		for (size_t i = l - 1 ; i >= 1 ; --i)
		{	u.push_back(*new grid (i));
			r.push_back(*new grid (i));
			e.push_back(*new grid (i));
		}
	}
	
	~mgsolve()
	{}

	void solve()
	{	double current_qs = 1;
		double last_qs = 1;
		double norm_residual = 0;
		double norm_error = 0;
		double conv_rate = 0;

		e.begin()->zero();
		
		struct timeval tim;
		gettimeofday(&tim , NULL);
		double t_start = tim.tv_sec * 1000.0 + (tim.tv_usec / 1000.0);
		
		compute_residual (e.begin() , u.begin() , r.begin());
		current_qs = r.begin()->quad_sum();

		for (size_t i = 0 ; i < n ; ++i)
		{	mgm (u.begin() , e.begin());
			
			compute_residual (e.begin() , u.begin() , r.begin());
			
			last_qs = current_qs;
			current_qs = r.begin()->quad_sum();
			
			conv_rate = sqrt(current_qs / last_qs);	//convergence rate
			norm_residual = sqrt (current_qs) / (r.begin()->get_width() * r.begin()->get_height());	//L2 norm of residual
			
			cout << "V-cycle " << i << ":" << endl;
			cout << "L2 norm of residual is: " << norm_residual << ";" << endl;
			cout << "Convergence rate is: " << conv_rate << "." << endl << endl;
		}
		
		gettimeofday(&tim , NULL);
		double t_end = tim.tv_sec * 1000.0 + (tim.tv_usec / 1000.0);
		cout << "Time used: " << t_end - t_start << " ms." << endl << endl;
		
		grid e = diff(exact_solution() , *u.begin());
		norm_error = sqrt (e.quad_sum()) / (e.get_width() * e.get_height());
		cout << "L2 norm of error is: " << norm_error << "." << endl << endl;

		cout << "Outputing final approximation into file \"solutin.txt\" ..." << endl;
		u.begin() -> output();
		cout << "Output completed!" << endl << endl;
	}

	grid exact_solution()
	{	double h = 1 / pow (2 , l);
		grid u_star (l);

		for (long i = 0 ; i < pow (2 , l) + 1 ; ++i)
		{	for (long j = 0 ; j < pow (2 , l) + 1 ; ++j)
			{	double temp = sin(M_PI * h * i) * sinh (M_PI * h * j);
				u_star.set (i , j , temp);
			}
		}
		return u_star;
	}
};



int main (int argc, char* argv[]){
	if(argc != 3){
	   cerr << "Illegual call" << endl;
	   cerr << "Required: ./mgsolve level num_of_cycles " << endl;
	   return -1;
	}
	
	int l = atoi(argv[1]);
	int n = atoi(argv[2]);
	
	mgsolve test(l , n);
	test.solve();
	
	return 0;
}




