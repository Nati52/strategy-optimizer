#include <cxxabi.h>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

using namespace std;

typedef vector<int> vi;
typedef pair<int, int> pii;
typedef pair<double, double> pdd;

#define eat(a) int a; cin >> a
#define eatd(a) double a; cin >> a
#define eb emplace_back
#define f first
#define s second
#define f0r(i, a) for (int i = 0; i < a; i++)
#define f1r(i, a, b) for (int i = a; i < b; i++)

struct depend {
	vector<int> reqs;
	vector<pii> crit;
	vector<int> most;
	int sum;
	int val;
};

const double comp_time = 150;
const double shop_hour = 31 * 6;

int n, n2, m;
int ma = 0;

map<char, int> dict;			// maps chars to mechanisms
vector<vector<char> > mechs;	// list of mechs for each action/arch_type
vi rep;							// value of inf rep
vector<vi> total;				// list of total points for each action (before inf rep)
vi len;						// max times you can repeat, or 1e7 if there's inf rep
vector<pdd> t;				// effective time for each action
vector<double> arch_t; 		// archtype times (may need to input t, q instead of just t)

vector<vector<char> > mechs2;	// list of mechs for each action/arch_type
vi rep2;						// value of inf rep
vector<vi> total2;				// list of total points for each action (before inf rep)
vi len2;						// max times you can repeat, or 1e7 if there's inf rep
vector<pdd> t2;					// effective time for each action
vector<double> arch_t2; 		// archtype times (may need to input t, q instead of just t)

vector<double> mech_bt; 	// mechanism build times

vi freq;			// how often each action is performed
vi freqy;			// how often each criteria is performed
vector<bool> in; 	// whether we've used mechanism yet
vi upper;

vector<vi> best;	// optimal distribution
vector<vi> best2; 	// optimal distribution
vector<vector<pii> > pairs;
vi score; 			// our score

vector<depend> deps; // Dependency info

int point(int i) {
	if (freq[i] >= total[i].size()) {
		int last = *(total[i].end() - 1);
		return last + (freq[i] - total[i].size() + 1) * rep[i];
	}
	return total[i][freq[i]];
}

void iterate(int ind, double et_left, double bt_left, double at1_left, double at2_left) {		// Assumes last action can be repeated infinitely, this is easy to do in the input
	double delta_bt = 0.0;
	vector<char> to_add;
	for (char c : mechs[ind]) {
		if (!in[dict[c]]) {
			to_add.eb(c);
			delta_bt += mech_bt[dict[c]];
		}
	}

	int our_max = min(1.0 * len[ind], et_left / t[ind].f);
	if (delta_bt > bt_left)
		our_max = 0;

	if (ind == n - 1) {
		freq.eb(our_max);

		int ans = 0;
		f0r (i, n)
			ans += point(i);

		f0r (i, n2)
			ans += rep2[i] * freqy[i];

		int us = ans;
		ans *= n * (n + 1) / 2;

		vector<pii> arch_pairs;
		int temp_max = 0;

		f0r (i, n) f1r(j, i, n) {
				int change = 0;
				if (i == j) {
					int p = min(1.0 * len[i] - freq[i], floor(at1_left / arch_t[i]) + floor(at2_left / arch_t[i]));
					freq[i] += p;
					change += point(i);
					for (depend d : deps) {
						int sum = 0;
						for (int req : d.reqs)
							sum += freq[req];
						if (sum < d.sum)
							continue;

						bool ok = 1;
						for (pii crit : d.crit)
							if (freqy[crit.f] != crit.s)
								ok = false;

						if (ok)
							ans += d.val;
					}
					freq[i] -= p;
					change -= point(i);
				}
				else {
					int p1 = min(1.0 * len[i] - freq[i], floor(at1_left / arch_t[i]));
					int p2 = min(1.0 * len[j] - freq[j], floor(at2_left / arch_t[j]));
					freq[i] += p1;
					freq[j] += p2;
					change += point(i);
					change += point(j);
					for (depend d : deps) {
						int sum = 0;
						for (int req : d.reqs)
							sum += freq[req];
						if (sum < d.sum)
							continue;

						bool ok = 1;
						for (pii crit : d.crit)
							if (freqy[crit.f] != crit.s)
								ok = false;

						if (ok)
							ans += d.val;
					}
					freq[i] -= p1;
					freq[j] -= p2;
					change -= point(i);
					change -= point(j);
				}
				if (change > temp_max) {
					temp_max = change;
					arch_pairs.resize(0);
				}
				if (change >= temp_max)
					arch_pairs.eb(i, j);
				ans += change;
			}

		if (ans > ma) {
			best.resize(0);
			best2.resize(0);
			score.resize(0);
			pairs.resize(0);
			ma = ans;
		}
		if (ans >= ma) {
			best.eb(freq);
			best2.eb(freqy);
			score.eb(us);
			pairs.eb(arch_pairs);
		}

		freq.erase(freq.end() - 1);
		return;
	}
	else {
		f0r (i, our_max + 1) {
			freq.eb(i);
			if (i == 0) {
				iterate(ind + 1, et_left, bt_left, at1_left, at2_left);
				for (char c : to_add)
					in[c] = true;
			}
			else
				iterate(ind + 1, et_left - i * t[ind].f, bt_left - delta_bt, at1_left, at2_left);
			freq.erase(freq.end() - 1);
		}

		for (char c : to_add)
			in[c] = false;
	}
}

void iterate2(int ind, double et_left, double bt_left, double at1_left, double at2_left) {
	f0r (i, upper[ind] + 1) {
		double delta_bt = 0.0;
		vector<char> to_add;
		if (i) {
			for (char c : mechs2[ind]) {
				if (!in[dict[c]]) {
					to_add.eb(c);
					delta_bt += mech_bt[dict[c]];
				}
			}
			bt_left -= delta_bt;
			for (char c : to_add)
				in[c] = true;
		}

		f0r (j, upper[ind] + 1) f0r (k, upper[ind] + 1) {
				freqy.eb(i + j + k);
				if (ind != n2 - 1)
					iterate2(ind + 1, et_left - t2[ind].f * i, bt_left, at1_left - arch_t2[ind] * j, at2_left - arch_t2[ind] * k);
				else
					iterate(0, et_left - t2[ind].f * i, bt_left, at1_left - arch_t2[ind] * j, at2_left - arch_t2[ind] * k);
				freqy.erase(freqy.end() - 1);
			}
		if (i)
			for (char c : to_add)
				in[c] = false;
	}
}

int main() {
	cin >> n >> n2 >> m;

	f0r (i, m) {
		char c;
		scanf(" %c ", &c);
		dict[c] = i;
	}
	f0r (i, m) {
		eat(bt);
		mech_bt.eb(bt);
	}

	rep.resize(n);
	total.resize(n);
	mechs.resize(n);
	in.resize(n);

	rep2.resize(n2);
	total2.resize(n2);
	mechs2.resize(n2);
	upper.resize(n2);

	f0r (i, n) {
		total[i].eb(0);
		eat(num_its);
		bool inf = false;

		f0r (j, num_its) {
			eat(val);

			if (val == -1)
				inf = true;

			int last = total[i][total[i].size() - 1];
			if (val < 0) {
				eat(rep_points);
				j++;

				if (val == -1) 		// Infinite repetition
					rep[i] = rep_points;
				else { 				// Finite repetition
					int num_reps = -val;
					f0r (k, num_reps)
						total[i].eb(last + (k + 1) * rep_points);
				}
			}
			else 					// No repetition
				total[i].eb(last + val);
		}

		if (inf)
			len.eb(1e7);
		else
			len.eb(total[i].size() - 1);
	}

	f0r (i, n) {
		eat(num_mechs);
		double et_eff = 0.0;
		double bt_eff = 0.0;

		f0r (j, num_mechs) {
			eatd(exec_time);
			char c;
			scanf(" %c ", &c);
			eatd(competence);

			mechs[i].eb(c);
			et_eff += exec_time / competence;	// Assumes will not have to restart after one mechanism fails.
			bt_eff += mech_bt[dict[c]]; 		// If not, then we'll have to do some math...
		}

		t.eb(et_eff, bt_eff);
	}

	f0r (i, n) {
		eatd(t_i);
		arch_t.eb(t_i);
	}


	f0r (i, n2) {
		total2[i].eb(0);
		eat(num_its);
		bool inf = false;

		f0r (j, num_its) {
			eat(val);

			if (val == -1)
				inf = true;

			int last = total2[i][total2[i].size() - 1];
			if (val < 0) {
				eat(rep_points);
				j++;

				if (val == -1)		// Infinite repetition
					rep2[i] = rep_points;
				else { 				// Finite repetition
					int num_reps = -val;
					f0r (k, num_reps)
						total2[i].eb(last + (k + 1) * rep_points);
				}
			}
			else 					// No repetition
				total2[i].eb(last + val);
		}

		if (inf)
			len2.eb(1e7);
		else
			len2.eb(total2[i].size() - 1);
	}

	f0r (i, n2) {
		eat(num_mechs);
		double et_eff = 0.0;
		double bt_eff = 0.0;

		f0r (j, num_mechs) {
			eatd(exec_time);
			char c;
			scanf(" %c ", &c);
			eatd(competence);

			mechs2[i].eb(c);
			et_eff += exec_time / competence;	// Assumes will not have to restart after one mechanism fails.
			bt_eff += mech_bt[dict[c]]; 		// If not, then we'll have to do some math...
		}

		t2.eb(et_eff, bt_eff);
	}

	f0r (i, n2) {
		eatd(t_i);
		arch_t2.eb(t_i);
	}

	// Add in number of times criteria needs to be completed to input
	eat(num_deps);
	f0r (i, num_deps) {
		depend d;

		eat(sum);
		eat(val);
		d.sum = sum;
		d.val = val;

		eat(num_reqs);
		f0r (j, num_reqs) {
			eat(act);
			d.reqs.eb(act);
		}
		eat(num_crit);
		f0r (j, num_crit) {
			eat(atmost);
			eat(act);
			eat(num);
			d.crit.eb(act, num);
			d.most.eb(atmost);
			upper[act] = atmost;
		}

		deps.eb(d);
	}

	// Output the cummulative point values
	f0r (i, n) {
		cout << t[i].f << ' ' << t[i].s << endl;
		f0r (j, total[i].size())
			cout << total[i][j] << ' ';
		cout << endl;
	}
	cout << "--------------" << endl;

	iterate2(0, comp_time, shop_hour, comp_time, comp_time);

	cout << "Best Average Alliance Score:  " << ma * 2 / n / (n + 1) << endl;

	f0r (num, best.size()) {
		cout << "Your score: " << score[num] << endl << "Frequency List: ";
		f0r (i, n)
			cout << best[num][i] << ' ';
		cout << endl << "Freqy List: ";
		f0r (i, n2)
			cout << best2[num][i] << ' ';
		cout << endl << "Best pairs: ";
		for (pii x : pairs[num])
			cout << x.f + 1 << ' ' << x.s + 1 << '\t';
		cout << endl;
	}
}