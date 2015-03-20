#include <set>
#include <tuple>

using namespace std;

class Retta 
{
private:
  tuple<int,int> p[3];
  int prob;

  void calcola_prob()
  {
    /* TODO */
    this->prob = 1;
  }
public:
  Retta(tuple<int,int> *inp)
  {
    for (int i=0; i<3; ++i)
      {
			this->p[i] = inp[i];
      }
    calcola_prob();
  }
  inline int get_prob()
  {
    return prob;
  }
  inline get_points(tuple<int,int> outp)
  {
    for (int i=0; i<3; ++i)
      {
	outp[i] = this->p[i];
      }
  }
  friend bool operator< (const Retta& r, const Retta& s)
  {
    if (r.prob != s.prob)
      {
	return (r.prob < s.prob);
      }
    else
      {
	// TODO: per maggior efficienza andrebbero confrontati le triplette
	return true;
      }
  }
};

void gioca()
{
  // input
  vector<tuple<int,int>> punti[3];

  // data
  set<Retta> tutto; //tutti le n^3 rette
  map<tuple<int,int>,set<set<Retta>::iterator> > riferimenti[3]; // per ricordarsi dove compare un certo punto
  for (vector<tuple<int,int> >::iterator p0=punti[0].begin(); p0 != punti[0].end() ; ++i)
    {
      for (vector<tuple<int,int> >::iterator p1=punti[0].begin(); p0 != punti[0].end() ; ++i)
	{
	  for (vector<tuple<int,int> >::iterator p2=punti[0].begin(); p0 != punti[0].end() ; ++i)
	    {
	      tuple<int,int> inp[3] = { *p0, *p1, *p2 };
	      // Inserisco la retta nel set
	      set<Retta>::iterator pos = tutto.insert(Retta(inp));
	      if (pos->get_prob()>0)
		{
		  // e` un punto possibile, allora metto i riferimenti
		  for (int i=0; i<3; ++i)
		    {
		      riferimenti[i][inp[i]].insert(pos);
		    }
		}
	      else
		{
		  // e` un evento impossibile, lo butto via
		  tutto.erase(pos);
		}
	    }
	}
    }
  vector<Retta> uscita; // i punti selezionati
  while (!tutto.empty())
    {
      set<Retta>::iterator pos = tutto.upper_bound();
      uscita.push_back(*pos);
      tuple<int,int> p[3];
      pos->get_points(p);
      // elimina le rette che usano questi punti
      for (int i=0; i<3; ++i)
	{
	  for (set<set<Retta>::iterator>::iterator j=riferimenti[i][p[i]].begin(); j!=riferimenti[i][p[i]].end(); ++j)
	    {
	      tutto.erase(*j);
	    }
	}
    }
}
