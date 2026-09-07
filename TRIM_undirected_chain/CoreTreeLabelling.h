#ifndef CORETREELABELLING_H_
#define CORETREELABELLING_H_

#include <omp.h>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <queue>
#include <bitset>
#include <cstring>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <climits>
#include <iostream>
#include <functional>
#include <chrono>
#include <iomanip>
#include <sys/time.h>

using namespace std;

#define MAXN 1000000000
#define MAXINT ((unsigned)4294967295)
#define MAXGID 100000000
#define MAXLINE 1024
#define MAXT 65000
#define MAXD 120
#define POSINF 2147483647
#define NEGINF 0

#define N_ROOTS 4
#define MAX_BP_THREADS 8

#define RANK_STATIC 0
#define RANK_LOCAL_STATIC 1	   // LS
#define RANK_HOP_BETWEENNESS 2 // HB

typedef unsigned short tint;
typedef char dint;

vector<int> score, _score;

class VAL
{
public:
	VAL(int id, int val);
	int id, val;
	bool operator<(const VAL &v) const;
};
bool VAL::operator<(const VAL &v) const
{
	return val == v.val ? id < v.id : val > v.val;
}

VAL::VAL(int id, int val)
{
	this->id = id;
	this->val = val;
}
double get_time(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1e-6;
}
class DIS
{
public:
	DIS() : id(0), dis() {} // ??????
	DIS(int id, vector<int> dis);
	int id;
	vector<int> dis;
	bool operator<(const DIS &v) const;
	void clear()
	{
		dis.clear(); // ??dis
		id = 0;		 // ?????????id??????
	}
};

DIS::DIS(int id, vector<int> dis) : id(id), dis(dis) {}

bool DIS::operator<(const DIS &v) const
{
	return id < v.id;
}

class IntVal
{
public:
	int x;

public:
	IntVal();
	IntVal(int _x);
	bool operator<(const IntVal &v) const;
};

class DV
{
public:
	DV();
	DV(int id, double val);
	int id;
	double val;
	bool operator<(const DV &v) const;
};

class CM
{
public:
	CM() : id(0) {}
	CM(int id, int fir);
	int id;
	vector<int> ve;
};

class Vertex
{
public:
	vector<int> adj;
	int lb, ub, len;
};

class Val
{
public:
	Val(int id, int val);
	int id, val;
	bool operator<(const Val &v) const;
};
bool Val::operator<(const Val &v) const
{
	return val == v.val ? id < v.id : val > v.val;
}
Val::Val(int id, int val)
{
	this->id = id;
	this->val = val;
}

class MEM
{
public:
	int u, source, dis;

	MEM() : u(0), source(0) {}

	MEM(int _u, int _source) : u(_u), source(_source) {}
};
struct MEM1
{
	int cur, source, dis;
	bool operator<(const MEM1 &other) const
	{
		return dis > other.dis;
	}
};
class TreeNode
{
public:
	int id, f, h, rid, rsize, w, ch_n;
	vector<int> nbr, cost; // of size <=w
	vector<int> anc;	   // of size h
	vector<dint> dis;
	vector<pair<int, int>> dpairs;
	vector<int> ch;
};

struct BPLabel
{
	uint8_t bpspt_d[N_ROOTS];
	uint64_t bpspt_s[N_ROOTS][2];
};

class CoreTree
{
public:
	static inline bool get_edge(char *line, int &a, int &b, int num_cnt = 2);
	static inline int get_num_cnt(string path);
	static void create_bin(string path, string graph_name, int rank_threads = 1,
						   int rank_method = RANK_STATIC, int rank_max_minutes = 1000000, int max_hops = 3, bool merge_equv = true, string bin_path = "");
	static void get_order(vector<int> *con, int n, int *o, int method);
	static void get_order_ls(vector<int> *con, int n, int *o, int rank_threads, int rank_max_minutes);
	static void get_order_hb(vector<int> *con, int n, int *o, int rank_threads, int rank_max_minutes, int n_hop);

public:
	string path;
	int **con, *dat, *deg, *nid, *old, *old_new, *new_old, *f;
	int n_org, n, n_core;
	long long m, m_core;
	unsigned MAXDIS, MAXMOV, MASK;
	int root_n;

public:
	CoreTree(string path);
	~CoreTree();
	void load_graph(string bin_path);

public:
	int id;
	double t;
	int max_w;
	int numLand;
	int stamp_vis, stamp_dis, stamp, nnn;
	vector<int> d_l;
	vector<vector<int>> nbr, cost;
	vector<int> ord, rank;
	vector<int> root;
	vector<int> root_to_idex;
	vector<vector<pair<int, int>>> E;
	vector<TreeNode> tree;
	vector<Vertex> V;
	vector<vector<int>> nodeGroup;
	vector<int> landId;
	vector<vector<int>> disLand;
	vector<vector<Val>> listLand;
	vector<int> father;
	vector<vector<DIS>> root_dis;
	vector<int> candidate;
	vector<vector<vector<int>>> sum_update_dis;
	vector<int> vi;
	vector<int> id_to_index;
	vector<int> stamp_d;
	vector<int> stamp_v;
	vector<int> Dijsktra_dis;
	vector<vector<int>> Muti_Dijsktra_dis;
	vector<vector<int>> Muti_Dijsktra_vis;
	vector<vector<int>> update_dis;
	vector<vector<int>> BFS_dis;
	vector<int> node_g;
	vector<int> group_l;
	vector<int> group_r;
	vector<int> group_len;
	vector<int> group_size;
	vector<vector<int>> group_d;
	vector<vector<int>> group_mem;
	vector<vector<int>> Dij_dis;
	vector<int> record_dis;
	vector<int> used_dis;
	vector<int> v_record;
	vector<int> v_used;

public:
	void reduce(int max_w, int n_threads);
	void create_tree();

	void decompose_tree(int max_w, int n_threads, string bin_path);

	void load_label_tree(int max_w, string bin_path);
	void load_label(int max_w, string bin_path, int numLand);

	void save_label_tree(int max_w, string bin_path);

	void init(int max_w);
	void create_landmark();
	void select_landmark();
	void distance_landmark();
	inline void update_ub(int v, int ub);
	inline void update_lb(int v, int lb);
	void update_bound();
	void assign_group();
	int get_unmch(int lid);
	void compute_ecc();
	void compute_ecc(int lid);
	void txt_to_file(string graph_name);

	void save_tmp_graph(int max_w, string bin_path);
	void load_tmp_graph(int max_w, string bin_path);
	void Compute_group();

public:
	int *dis;

public:
	void BFS(int u, vector<int> &distance);
	bool BFS_Judge(int u, vector<int> &node_vis, vector<int> &node_continue);
	void Local_Djikstra(int u, int stamp);
	void Local_Djikstra(int u, int stamp, vector<vector<int>> &up_dis);
	void Muti_Dijsktra(int u);
	void Local_BFS(int u, int lambda, int lid);
	void Compute_root_dis(int u);
	void Merge_root_dis();
	void update(int u, int v);
	bool check(int u, vector<int> &visited);
	bool short_dij(int x, int y, int dis_xy);
	void Count_branch_dis(int u, int v, vector<int> &dis, int &cnt);
	void init_pair(int u, int i);
	int update_pair(int u, int i);
	int binarySearch(const vector<int> &arr, int target, int g_len, int d1, int d2);
};

bool DV::operator<(const DV &v) const
{
	if (val > v.val + 1e-8)
		return true;
	if (val < v.val - 1e-8)
		return false;
	return id < v.id;
}

DV::DV(int id, double val)
{
	this->id = id;
	this->val = val;
}

struct node
{
	int u, d;
	node() {}
	node(int u, int d) : u(u), d(d) {}
	bool operator<(const node other) const
	{
		return d > other.d;
	}
};

DV::DV()
{
	id = -1;
	val = -1;
}

IntVal::IntVal() { x = -1; }
IntVal::IntVal(int x) { this->x = x; }
bool IntVal::operator<(const IntVal &v) const
{
	if (score[x] == score[v.x])
		return x < v.x;
	return score[x] < score[v.x];
}

CoreTree::CoreTree(string path)
{
	deg = NULL;
	dat = NULL;
	con = NULL;
	nid = NULL;
	old = NULL;
	f = NULL;
	n_org = 0;
	n = 0;
	m = 0;
	n_core = 0;
	m_core = 0;
	MAXDIS = 0;
	MAXMOV = 0;
	MASK = 0;
	max_w = 0;
	dis = NULL;
	this->path = path;
	t = omp_get_wtime();
}

CoreTree::~CoreTree()
{
	if (deg)
		delete[] deg;
	if (dat)
		delete[] dat;
	if (con)
		delete[] con;
	if (nid)
		delete[] nid;
	if (old)
		delete[] old;
	if (f)
		delete[] f;
	if (dis)
		delete[] dis;
}

void CoreTree::BFS(int u, vector<int> &distance)
{
	queue<int> toVisit;
	vector<int> visited;
	visited.resize(n, 0);

	int type = 0;
	int m = 0;

	toVisit.push(u);
	distance[u] = 0;
	visited[u] = 1;

	while (!toVisit.empty())
	{
		int cur = toVisit.front();
		toVisit.pop();
		m++;
		for (int i = 0; i < deg[cur]; i++)
		{
			int v = con[cur][i];
			if (visited[v] == 0)
			{
				visited[v] = 1;
				distance[v] = distance[cur] + 1;
				toVisit.push(v);
			}
		}
	}
}

bool CoreTree::BFS_Judge(int u, vector<int> &node_vis, vector<int> &node_continue)
{
	queue<MEM1> pq;
	vector<int> nbr;
	bool flag_chain;
	int group_index = node_g[u], cnt = 0, Dij_size = 0;
	int left = group_l[group_index], right = group_r[group_index];
	stamp_vis++;
	stamp_dis++;
	record_dis[group_mem[group_index].size()]++;
	if (deg[u] == 2 || deg[u] == 1)
	{
		if (node_continue[left] == 1 && node_continue[right] == 1)
			return true;
		if (!node_continue[left])
		{
			pq.push(MEM1{left, Dij_size, 0});
			Dij_dis[left].resize(n);
			Dij_dis[left][left] = 0;
			Muti_Dijsktra_dis[Dij_size][left] = stamp_dis;
			nbr.push_back(left);
			Dij_size++;
		}

		if (!node_continue[right] && left != right)
		{
			pq.push(MEM1{right, Dij_size, 0});
			Dij_dis[right].resize(n);
			Dij_dis[right][right] = 0;
			Muti_Dijsktra_dis[Dij_size][right] = stamp_dis;
			nbr.push_back(right);
			Dij_size++;
		}
		flag_chain = true;
		left = group_mem[group_index].front(), right = group_mem[group_index].back();
	}
	else
	{

		pq.push(MEM1{u, 0, 0});
		Dij_dis[u].resize(n);
		Dij_dis[u][u] = 0;
		Muti_Dijsktra_dis[Dij_size][u] = stamp_dis;
		nbr.push_back(u);
		flag_chain = false;
		Dij_size++;
	}
	while (!pq.empty())
	{
		MEM1 member = pq.front();
		int cur = member.cur;
		int source = member.source;
		int nowDis = member.dis;
		pq.pop();
		if (Muti_Dijsktra_vis[source][cur] == stamp_vis)
		{
			continue;
		}
		Muti_Dijsktra_vis[source][cur] = stamp_vis;
		if (flag_chain == false)
		{
			if (node_continue[cur] && d_l[cur] >= nowDis)
			{
				{
					node_vis[u] = 1;
					if (!node_continue[u])
						vector<int>().swap(Dij_dis[u]);
					return false;
				}
			}
		}
		else
		{
			for (int i = 0; i < Dij_size; i++)
			{
				if (source == i || (Muti_Dijsktra_vis[i][cur] != stamp_vis && !node_continue[nbr[i]] && left != right))
					continue;
				if (node_continue[cur] && (d_l[cur] >= Dij_dis[nbr[0]][cur] || d_l[cur] >= Dij_dis[nbr[1]][cur]))
				{
					int u_index = binarySearch(group_d[group_index], (group_len[group_index] + Dij_dis[right][cur] - Dij_dis[left][cur]) / 2, group_len[group_index], Dij_dis[nbr[0]][cur], Dij_dis[nbr[1]][cur]);
					int d = min(Dij_dis[left][cur] + group_d[group_index][u_index], Dij_dis[right][cur] + group_len[group_index] - group_d[group_index][u_index]);
					if (d_l[cur] >= d)
					{
						for (int i = 0; i < group_mem[group_index].size(); i++)
						{
							node_vis[group_mem[group_index][i]] = 1;
						}
						if (!node_continue[left])
							vector<int>().swap(Dij_dis[left]);
						if (!node_continue[right])
							vector<int>().swap(Dij_dis[right]);
						return false;
					}
				}
			}
		}
		for (int i = 0; i < deg[cur]; i++)
		{
			int d = nowDis + 1;
			int v = con[cur][i];
			if (Muti_Dijsktra_dis[source][v] != stamp_dis)
			{
				Muti_Dijsktra_dis[source][v] = stamp_dis;
				Dij_dis[nbr[source]][v] = d;
				pq.push(MEM1{v, source, d});
			}
		}
	}
	if (deg[u] == 2)
	{
		for (int i = 0; i < group_mem[group_index].size(); i++)
		{
			node_vis[group_mem[group_index][i]] = 1;
		}
		node_continue[left] = node_continue[right] = 1;
	}
	else
	{
		node_vis[u] = 1;
		node_continue[u] = 1;
	}
	return true;
}
void CoreTree::Local_Djikstra(int u, int stamp)
{

	int nbr_size = tree[tree[u].rid].nbr.size();
	vector<queue<int>> bucket;

	int cnt = 0;
	bucket.resize(2);
	bucket[0].push(u);
	Dijsktra_dis[u] = 0;
	stamp_dis++;
	stamp_vis++;
	stamp_d[u] = stamp_dis;

	while (cnt < bucket.size())
	{
		if (bucket[cnt].empty())
		{
			cnt++;
			continue;
		}
		int cur = bucket[cnt].front();
		bucket[cnt].pop();
		if (stamp_v[cur] == stamp_vis)
			continue;
		stamp_v[cur] = stamp_vis;
		if (vi[cur] == stamp && nbr_size != 1)
		{
			for (int i = 0; i < nbr_size; i++)
			{
				int v = tree[tree[u].rid].nbr[i];
				if ((stamp_d[v] != stamp_dis) || (Dijsktra_dis[v] > Dijsktra_dis[cur] + update_dis[id_to_index[cur]][tree[tree[u].rid].nbr[i]]))
				{
					Dijsktra_dis[v] = Dijsktra_dis[cur] + update_dis[id_to_index[cur]][tree[tree[u].rid].nbr[i]];
					stamp_d[v] = stamp_dis;
					if (Dijsktra_dis[v] >= bucket.size())
					{
						bucket.resize(Dijsktra_dis[v] + 1);
					}
					bucket[Dijsktra_dis[v]].push(v);
				}
			}
		}
		for (int i = 0; i < deg[cur]; i++)
		{
			int v = con[cur][i];
			int d = 1;
			if ((stamp_d[v] != stamp_dis || Dijsktra_dis[v] > Dijsktra_dis[cur] + d) && (rank[v] != -1 || vi[v] == stamp))
			{
				Dijsktra_dis[v] = Dijsktra_dis[cur] + d;
				stamp_d[v] = stamp_dis;
				if (Dijsktra_dis[v] >= bucket.size())
				{
					bucket.resize(Dijsktra_dis[v] + 1);
				}
				bucket[Dijsktra_dis[v]].push(v);
			}
		}
	}
}

void CoreTree::Local_Djikstra(int u, int stamp, vector<vector<int>> &up_dis)
{

	int nbr_size = tree[tree[u].rid].nbr.size();
	vector<queue<int>> bucket;

	int cnt = 0;
	bucket.resize(2);
	bucket[0].push(u);
	Dijsktra_dis[u] = 0;
	stamp_dis++;
	stamp_vis++;
	stamp_d[u] = stamp_dis;

	while (cnt < bucket.size())
	{
		if (bucket[cnt].empty())
		{
			cnt++;
			continue;
		}
		int cur = bucket[cnt].front();
		bucket[cnt].pop();

		if (stamp_v[cur] == stamp_vis)
			continue;
		stamp_v[cur] = stamp_vis;
		if (vi[cur] == stamp && nbr_size != 1)
		{
			for (int i = 0; i < nbr_size; i++)
			{
				int v = tree[tree[u].rid].nbr[i];
				if ((stamp_d[v] != stamp_dis) || (Dijsktra_dis[v] > Dijsktra_dis[cur] + up_dis[id_to_index[cur]][i]))
				{
					Dijsktra_dis[v] = Dijsktra_dis[cur] + up_dis[id_to_index[cur]][i];
					stamp_d[v] = stamp_dis;
					if (Dijsktra_dis[v] >= bucket.size())
					{
						bucket.resize(Dijsktra_dis[v] + 1);
					}
					bucket[Dijsktra_dis[v]].push(v);
				}
			}
		}
		for (int i = 0; i < deg[cur]; i++)
		{
			int v = con[cur][i];
			if (!(rank[v] != -1 || vi[v] == stamp))
				continue;

			int d = 1;
			if (stamp_d[v] != stamp_dis || Dijsktra_dis[v] > Dijsktra_dis[cur] + d)
			{
				Dijsktra_dis[v] = Dijsktra_dis[cur] + d;
				stamp_d[v] = stamp_dis;
				if (Dijsktra_dis[v] >= bucket.size())
				{
					bucket.resize(Dijsktra_dis[v] + 1);
				}
				bucket[Dijsktra_dis[v]].push(v);
			}
		}
	}
}
void CoreTree::Muti_Dijsktra(int u)
{
	vector<queue<MEM>> bucket;
	vector<int> nbr;
	vector<int> sum_result;
	vector<vector<int>> nbr_vis;

	int nbr_size = tree[tree[u].rid].nbr.size();

	if (nbr_size == 1)
		return;
	int cnt = 0, sum = 0;
	nbr.assign(nbr_size, 0);
	sum_result.resize(nbr_size);
	bucket.resize(2);
	nbr_vis.resize(nbr_size);
	stamp_vis++;
	stamp_dis++;

	for (int i = 0; i < nbr_size; i++)
	{
		nbr[i] = tree[tree[u].rid].nbr[i];
		nbr_vis[i].assign(nbr_size, 0);
	}

	for (int i = 0; i < nbr_size; i++)
	{
		update_dis[i][nbr[i]] = 0;
		Muti_Dijsktra_dis[i][nbr[i]] = stamp_dis;
		bucket[0].push(MEM{nbr[i], i});
	}

	if (bucket[cnt].size() == 0)
		return;

	while (1)
	{
		while (bucket[cnt].size() == 0)
			cnt++;

		MEM member = bucket[cnt].front();
		int cur = member.u;
		int source = member.source;
		bucket[cnt].pop();
		if (sum_result[source] == nbr_size)
			continue;

		if (Muti_Dijsktra_vis[source][cur] == stamp_vis)
		{
			continue;
		}
		Muti_Dijsktra_vis[source][cur] = stamp_vis;
		for (int i = 0; i < nbr_size; i++)
		{
			if (source == i || Muti_Dijsktra_vis[i][cur] != stamp_vis || nbr_vis[i][source] || sum_result[i] == nbr_size - 1)
				continue;
			update_dis[source][nbr[i]] = cnt + update_dis[i][cur];
			update_dis[i][nbr[source]] = cnt + update_dis[i][cur];
			sum_result[source]++;
			sum_result[i]++;
			nbr_vis[i][source] = nbr_vis[source][i] = 1;
		}

		for (int i = 0; i < deg[cur]; i++)
		{
			int d = cnt + 1;
			int v = con[cur][i];
			if ((Muti_Dijsktra_dis[source][v] != stamp_dis) || (update_dis[source][v] > d))
			{
				Muti_Dijsktra_dis[source][v] = stamp_dis;
				update_dis[source][v] = d;
				if (d >= bucket.size())
				{
					bucket.resize(d + 1);
				}
				bucket[d].push(MEM{v, source});
			}
		}

		sum = 0;
		for (int i = 0; i < nbr_size; i++)
		{
			sum += sum_result[i];
		}
		if (sum == nbr_size * (nbr_size - 1))
		{
			break;
		}
	}
}

void CoreTree::Local_BFS(int u, int lambda, int lid)
{
	queue<int> toVisit;

	int nbr_size = tree[tree[u].rid].nbr.size();
	int u_index = root_to_idex[tree[u].rid];
	int fa_index = father[u_index];

	stamp = ++stamp_dis;
	for (int i = 0; i < tree[tree[u].rid].nbr.size(); i++)
	{
		vi[tree[tree[u].rid].nbr[i]] = stamp;
		id_to_index[tree[tree[u].rid].nbr[i]] = i;
	}

	Local_Djikstra(u, stamp, sum_update_dis[fa_index]);

	for (int i = 0; i < root_dis[fa_index].size(); i++)
	{
		int v = root_dis[fa_index][i].id;
		int dis = Dijsktra_dis[v];

		int id = v;
		int lb = V[id].lb, ub = V[id].ub;
		if (lb < ub)
		{
			update_lb(id, dis);
			if (V[id].lb == V[id].ub)
				nnn--;
		}
	}
}

bool CoreTree::check(int u, vector<int> &visited)
{

	if (rank[u] == -1)
	{
		return false;
	}

	int u_index = root_to_idex[tree[u].rid];
	int fa_index = father[u_index];

	if (!visited[fa_index])
	{
		Compute_root_dis(u);
		visited[fa_index] = 1;
	}

	if (candidate[u])
	{
		return false;
	}
	else
		return true;
}

void CoreTree::Compute_root_dis(int u)
{
	queue<int> que;

	int v, d, father_index, cnt;
	father_index = father[root_to_idex[tree[u].rid]];
	for (int i = 0; i < root.size(); i++)
	{
		if (father[i] == father_index)
		{
			que.push(root[i]);
		}
	}
	Muti_Dijsktra(u);
	stamp = ++stamp_dis;
	for (int i = 0; i < tree[tree[u].rid].nbr.size(); i++)
	{
		vi[tree[tree[u].rid].nbr[i]] = stamp;
		id_to_index[tree[tree[u].rid].nbr[i]] = i;
	}

	// vector<vector<int>> tmp_dis;
	// int nbr_size = tree[tree[u].rid].nbr.size();
	// tmp_dis.resize(nbr_size);
	// for (int i = 0; i < nbr_size; i++)
	// {
	// 	for (int j = 0; j < nbr_size; j++)
	// 		tmp_dis[i].push_back(update_dis[i][tree[tree[u].rid].nbr[j]]);
	// }
	// sum_update_dis[father[root_to_idex[tree[u].rid]]] = tmp_dis;

	int nbr_size = tree[tree[u].rid].nbr.size();
	sum_update_dis[father[root_to_idex[tree[u].rid]]].resize(nbr_size);
	for (int i = 0; i < nbr_size; i++)
	{
		for (int j = 0; j < nbr_size; j++)
			sum_update_dis[father[root_to_idex[tree[u].rid]]][i].push_back(update_dis[i][tree[tree[u].rid].nbr[j]]);
	}

	while (!que.empty())
	{
		int cur = que.front();
		que.pop();
		DIS tmp;
		Local_Djikstra(cur, stamp);
		for (int j = 0; j < tree[tree[cur].rid].nbr.size(); j++)
		{
			v = tree[tree[cur].rid].nbr[j];
			d = Dijsktra_dis[v];
			tmp.dis.push_back(d);
		}
		tmp.id = cur;
		root_dis[father_index].push_back(tmp);
		for (int j = 0; j < tree[cur].ch.size(); j++)
		{
			v = tree[cur].ch[j];
			que.push(v);
		}
	}

	vector<bool> used;
	used.assign(root_dis[father_index].size(), 1);

	sort(root_dis[father_index].begin(), root_dis[father_index].end(), [](const DIS &lhs, const DIS &rhs)
		 { return lhs.dis > rhs.dis; });

	for (int i = 0; i < root_dis[father_index].size(); i++)
	{
		if (used[i] == 0)
			continue;
		candidate[root_dis[father_index][i].id] = 1;
		for (int j = i + 1; j < root_dis[father_index].size(); j++)
		{
			cnt = 0;
			if (used[j] == 0)
				continue;
			for (int k = 0; k < tree[root[father_index]].nbr.size(); k++)
			{
				if (root_dis[father_index][i].dis[k] < root_dis[father_index][j].dis[k])
					break;
				if (root_dis[father_index][i].dis[k] >= root_dis[father_index][j].dis[k])
					cnt++;
			}

			if (cnt == tree[root[father_index]].nbr.size())
				used[j] = 0;
		}
	}
}

void CoreTree::Merge_root_dis()
{
	father.resize(root_n);
	int p = 0, cnt = 0, nn = 0, flag;
	vector<CM> f;
	f.resize(root_n);
	for (int i = 0; i < root.size(); i++)
	{
		f[i].id = i;
		f[i].ve = tree[root[i]].nbr;
		sort(f[i].ve.begin(), f[i].ve.end());
	}

	sort(f.begin(), f.end(), [](const CM &a, const CM &b)
		 {
			if (a.ve.size() == b.ve.size()) {
				for (int i = 0; i < a.ve.size(); i++) {
					if (a.ve[i] == b.ve[i]) continue;
					return a.ve[i] > b.ve[i];
				}
			}
			return a.ve.size() > b.ve.size(); });

	if (f.size() != 0)
		father[f[0].id] = f[0].id;

	for (int i = 1; i < f.size(); i++)
	{
		flag = 0;
		if (f[p].ve.size() != f[i].ve.size())
		{
			father[f[i].id] = f[i].id;
			p = i;
		}
		else
		{
			for (int j = 0; j < f[p].ve.size(); j++)
			{
				if (f[p].ve[j] != f[i].ve[j])
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
			{
				father[f[i].id] = f[p].id;
				nn++;
			}
			else
			{
				father[f[i].id] = f[i].id;
				p = i;
			}
		}
	}
}

void CoreTree::save_tmp_graph(int max_w, string bin_path)
{
	// printf("Saving Tmp Graph...\n");

	char stw[16];
	sprintf(stw, "%d", max_w);

	FILE *fout = fopen((bin_path + "tmp-" + string(stw) + ".bin").c_str(), "wb");
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m_core, sizeof(long long), 1, fout);
	fwrite(rank.data(), sizeof(int), n, fout);
	fclose(fout);
	// printf("Tmp Graph Saved!\n");
}

void CoreTree::load_tmp_graph(int max_w, string bin_path)
{
	// printf("Loading Tmp Graph...\n");

	char stw[16];
	sprintf(stw, "%d", max_w);

	FILE *fin = fopen((bin_path + "tmp-" + string(stw) + ".bin").c_str(), "rb");
	fread(&n, sizeof(int), 1, fin);
	n_core = 0;
	fread(&m_core, sizeof(long long), 1, fin);
	rank.resize(n);
	fread(rank.data(), sizeof(int), n, fin);
	fclose(fin);
	// printf("Tmp Graph Loaded!\n");
}

void CoreTree::save_label_tree(int max_w, string bin_path)
{
	// printf("Saving Tree Label...\n");
	char stw[16];
	sprintf(stw, "%d", max_w);

	FILE *fout = fopen((bin_path + "label-tree-" + string(stw) + ".bin").c_str(), "wb");
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(rank.data(), sizeof(int), n, fout);

	fwrite(&root_n, sizeof(int), 1, fout);
	fwrite(root.data(), sizeof(int), root_n, fout);

	root_to_idex.resize(n, -1);
	for (int i = 0; i < root.size(); i++)
		root_to_idex[root[i]] = i;
	fwrite(root_to_idex.data(), sizeof(int), n, fout);

	for (int i = 0; i < n; ++i)
		if (rank[i] >= 0)
		{
			TreeNode &tn = tree[i];
			fwrite(&tn.rid, sizeof(int), 1, fout);
			fwrite(&tn.rsize, sizeof(int), 1, fout);
			fwrite(&tn.h, sizeof(int), 1, fout);
			fwrite(&tn.w, sizeof(int), 1, fout);
			fwrite(&tn.ch_n, sizeof(int), 1, fout);
			fwrite(tn.nbr.data(), sizeof(int), tn.w, fout);
			int ch_size = tn.ch.size();
			fwrite(tn.ch.data(), sizeof(int), tn.ch_n, fout);
		}
	fclose(fout);

	// printf("Tree Label Saved!\n");
}
int CoreTree::binarySearch(const vector<int> &arr, int target, int g_len, int d1, int d2)
{
	int u = -1;
	auto lower = upper_bound(arr.begin(), arr.end(), target);
	int maxx = 0;
	if (lower != arr.begin())
	{
		int diff = min(d1 + *(lower - 1), d2 + g_len - *(lower - 1));
		if (diff > maxx)
		{
			maxx = diff;
			u = (lower - 1) - arr.begin();
		}
	}
	if (lower != arr.end())
	{
		int diff = min(d1 + *lower, d2 + g_len - *lower);
		if (diff > maxx)
		{
			maxx = diff;
			u = lower - arr.begin();
		}
	}
	return u;
}
void CoreTree::Compute_group()
{
	vector<int> vis;
	int count = 0;
	vis.resize(n, 0);
	for (int i = 0; i < n; i++)
	{
		if (deg[i] == 1 || deg[i] > 2)
		{
			bool tmp_empty = false;
			node_g[i] = group_l.size();
			for (int j = 0; j < deg[i]; j++)
			{
				// cout << "i:" << i << ",j:" << j << endl;
				int right = con[i][j];
				if (vis[right] || deg[right] > 2 || deg[right] == 1)
					continue;
				vector<int> tmp = {i};
				vector<int> dd = {0};
				int d = 0, curr = right, prev = i;
				bool is_end = false;
				tmp_empty = true;
				d += 1;
				node_g[right] = group_l.size();
				vis[right] = 1;
				while (!is_end)
				{
					tmp.push_back(curr);
					dd.push_back(d);
					if (deg[curr] != 2)
					{
						is_end = true;
						group_l.push_back(i);
						group_r.push_back(curr);
						group_len.push_back(d);
						group_mem.push_back(tmp);
						group_d[count] = dd;
						count++;
					}
					else
					{
						for (int k = 0; k < deg[curr]; k++)
						{
							int next = con[curr][k];
							if (next != prev)
							{
								d += 1;
								prev = curr;
								curr = next;
								right = curr;
								node_g[right] = count;
								if (deg[right] == 2)
									vis[right] = 1;
								break;
							}
						}
					}
				}
			}
			if (!tmp_empty)
			{
				vector<int> tmp = {i};
				vector<int> dd = {0};
				group_l.push_back(i);
				group_r.push_back(i);
				group_len.push_back(1);
				group_d[count] = dd;
				group_mem.push_back(tmp);
				count++;
			}
		}
	}
	group_size.resize(group_l.size());
	for (int i = 0; i < group_l.size(); i++)
	{
		group_size[i] = group_mem[i].size();
	}
}

void CoreTree::load_label_tree(int max_w, string bin_path)
{
	// printf("Loading Tree Label...\n");
	dis = new int[n_core];

	char stw[16];
	sprintf(stw, "%d", max_w);

	FILE *fin = fopen((bin_path + "label-tree-" + string(stw) + ".bin").c_str(), "rb");

	fread(&n, sizeof(int), 1, fin);
	rank.resize(n);

	fread(rank.data(), sizeof(int), n, fin);
	n_core = 0;
	for (int i = 0; i < n; ++i)
		if (rank[i] == -1)
			++n_core;
	tree.resize(n);

	fread(&root_n, sizeof(int), 1, fin);
	root.resize(root_n);
	fread(root.data(), sizeof(int), root_n, fin);
	root_to_idex.resize(n);
	fread(root_to_idex.data(), sizeof(int), n, fin);
	int pos = 0;
	for (int i = 0; i < n; ++i)
		if (rank[i] >= 0)
		{
			TreeNode &tn = tree[i];
			fread(&tn.rid, sizeof(int), 1, fin);
			fread(&tn.rsize, sizeof(int), 1, fin);
			fread(&tn.h, sizeof(int), 1, fin);
			fread(&tn.w, sizeof(int), 1, fin);
			fread(&tn.ch_n, sizeof(int), 1, fin);
			tn.nbr.resize(tn.w);
			tn.anc.resize(tn.h);
			tn.dis.resize(tn.h);
			tn.ch.resize(tn.ch_n);
			fread(tn.nbr.data(), sizeof(int), tn.w, fin);
			// fread(tn.anc.data(), sizeof(int), tn.h-tn.w, fin);
			// fread(tn.dis.data(), sizeof(dint), tn.h, fin);
			fread(tn.ch.data(), sizeof(int), tn.ch_n, fin);
		}
	fclose(fin);
}

void CoreTree::load_label(int max_w, string bin_path, int numLand)
{
	this->max_w = max_w;
	this->numLand = numLand;
	load_label_tree(max_w, bin_path);
	load_tmp_graph(max_w, bin_path);
	load_graph(bin_path);
}

void CoreTree::decompose_tree(int max_w, int n_threads, string bin_path)
{
	if (con == NULL)
		load_graph(bin_path);
	double reduce_t = get_time();
	reduce(max_w, n_threads);
	create_tree();
	printf("%f\n", get_time() - reduce_t);
	save_label_tree(max_w, bin_path);
}

void CoreTree::init(int max_w)
{
	stamp_dis = 2;
	stamp_vis = 2;
	V.resize(n);
	root_dis.resize(root_n);
	candidate.resize(n, 0);
	sum_update_dis.resize(root.size());
	vi.resize(n);
	id_to_index.resize(n);
	stamp_d.assign(n, 0);
	stamp_v.assign(n, 0);
	Dijsktra_dis.assign(n, 0);
	Muti_Dijsktra_dis.resize(max_w);
	Muti_Dijsktra_vis.resize(max_w);
	update_dis.resize(max_w);
	d_l.resize(n, INT_MAX);
	node_g.resize(n, -1);
	group_d.resize(n);
	Dij_dis.resize(n);
	for (int i = 0; i < max_w; i++)
	{
		update_dis[i].assign(n, 0);
		Muti_Dijsktra_vis[i].assign(n, 0);
		Muti_Dijsktra_dis[i].assign(n, 0);
	}
	for (int i = 0; i < n; ++i)
	{
		sort(V[i].adj.begin(), V[i].adj.end());
		V[i].len = (int)V[i].adj.size();
		V[i].ub = POSINF;
		V[i].lb = NEGINF;
	}
}

void CoreTree::create_tree()
{
	tree.resize(n);
	for (int u = 0; u < n; ++u)
		tree[u].id = u;
	vector<pair<int, int>> v_pair;
	int maxh = 0, cnt_root = 0, maxdep = 0, max_sub_tree = 1;
	double tw = 0;
	for (int i = (int)ord.size() - 1; i >= 0; --i)
	{
		int x = ord[i];

		TreeNode &tn = tree[x];
		v_pair.clear();
		for (int j = 0; j < (int)nbr[x].size(); ++j)
		{
			int y = nbr[x][j];
			if (rank[y] == -1)
				v_pair.push_back(make_pair(n, j));
			else
				v_pair.push_back(make_pair(rank[y], j));
		}
		sort(v_pair.rbegin(), v_pair.rend());
		int w = (int)nbr[x].size();
		tn.nbr.resize(w);
		// tn.cost.resize(w);
		for (int j = 0; j < w; ++j)
		{
			tn.nbr[j] = nbr[x][v_pair[j].second];
			// tn.cost[j] = cost[x][v_pair[j].second];
		}

		tn.w = w;
		tn.id = x;
		tn.f = -1;
		tn.ch_n = 0;
		for (auto &u : nbr[x])
			if (rank[u] != -1 && (tn.f == -1 || rank[u] < rank[tn.f]))
				tn.f = u;
		if (tn.f == -1)
		{
			tn.h = tn.w + 1;
			++cnt_root;
			tn.rid = x;
			tn.rsize = tn.w;
			root.push_back(x);
		}
		else
		{
			tn.h = tree[tn.f].h + 1;
			tree[tn.f].ch.push_back(x);
			tree[tn.f].ch_n++;
			tn.rid = tree[tn.f].rid;
			tn.rsize = tree[tn.f].rsize;
		}
		tw += tn.rsize;
		maxh = max(maxh, tn.h);
	}
	root_n = cnt_root;

	// printf( "Core tree constructed, maxh=%d, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%0.3lf, t=%0.3lf secs\n",
	// maxh, maxdep, cnt_root, max_sub_tree, tw/(n-n_core), omp_get_wtime()-t);
}

void CoreTree::reduce(int max_w, int n_threads)
{
	// double t_t = omp_get_wtime();

	this->max_w = max_w;
	score.resize(n);
	_score.resize(n);
	set<IntVal> q;
	vector<vector<int>> bucket;
	vector<int> nodeToIndex;
	vector<int> nodeToDegree;
	nbr.resize(n);
	cost.resize(n);
	rank.resize(n);
	fill(rank.begin(), rank.end(), -1);
	bucket.resize(max_w);
	nodeToIndex.assign(n, -1);
	nodeToDegree.assign(n, -1);

	int r = 0;

	for (int i = 0; i < n; ++i)
		score[i] = _score[i] = deg[i];
	vector<bool> active(n, false);
	int cntt = 0;

	for (int u = 0; u < n; ++u)
	{
		if (deg[u] < max_w)
		{
			bucket[deg[u]].push_back(u);
			nodeToIndex[u] = bucket[deg[u]].size() - 1;
			nodeToDegree[u] = deg[u];
		}
		active[u] = true;
	}

	E.resize(n);

	for (int u = 0; u < n; ++u)
		for (int i = 0; i < deg[u]; ++i)
		{
			E[u].push_back(make_pair(con[u][i], 1));
		}

	int cnt = 0;

	vector<pair<int, int>> tmp;
	int nn = 0;
	while (nn < max_w)
	{
		if (nn >= max_w)
			break;
		if (bucket[nn].empty())
		{
			nn++;
			continue;
		}
		int x = bucket[nn].front();
		ord.push_back(x);
		int end_node = bucket[nn].back();
		if (end_node != x)
		{
			bucket[nodeToDegree[x]][nodeToIndex[x]] = end_node;
			nodeToDegree[end_node] = nodeToDegree[x];
			nodeToIndex[end_node] = nodeToIndex[x];
		}
		nodeToDegree[x] = -1;
		nodeToIndex[x] = -1;
		bucket[nn].pop_back();
		rank[x] = r++;
		for (auto &it : E[x])
		{
			nbr[x].push_back(it.first);
			cost[x].push_back(it.second);
		}

		for (auto &y : nbr[x])
		{
			if (E[y].size() >= max_w * 2)
			{
				active[y] = false;
			}
			if (!active[y])
				continue;
			for (size_t i = 0; i < E[y].size(); ++i)
			{
				if (E[y][i].first == x)
				{
					E[y].erase(E[y].begin() + i);
					break;
				}
			}
		}

		for (size_t i = 0; i < nbr[x].size(); ++i)
		{
			int u = nbr[x][i];
			if (!active[u])
			{
				E[u].push_back(make_pair(x, -cost[x][i]));
				continue;
			}
			tmp.clear();
			size_t j = 0, k = 0;
			while (j < nbr[x].size() && k < E[u].size())
				if (j == i)
					++j;
				else if (nbr[x][j] < E[u][k].first)
				{
					tmp.push_back(make_pair(nbr[x][j], cost[x][i] + cost[x][j]));
					++j;
				}
				else if (nbr[x][j] > E[u][k].first)
				{
					tmp.push_back(E[u][k]);
					++k;
				}
				else
				{
					if (E[u][k].second < cost[x][i] + cost[x][j])
						tmp.push_back(E[u][k]);
					else
						tmp.push_back(make_pair(nbr[x][j], cost[x][i] + cost[x][j]));
					++j;
					++k;
				}
			for (; j < nbr[x].size(); ++j)
				if (j != i)
					tmp.push_back(make_pair(nbr[x][j], cost[x][i] + cost[x][j]));
			for (; k < E[u].size(); ++k)
				tmp.push_back(E[u][k]);
			// update
			E[u] = tmp;
			_score[u] = (int)E[u].size();
			if (score[u] != _score[u])
			{
				score[u] = _score[u];
				if (nodeToIndex[u] != -1)
				{
					int end_node = bucket[nodeToDegree[u]].back();
					if (end_node != u)
					{
						bucket[nodeToDegree[u]][nodeToIndex[u]] = end_node;
						nodeToDegree[end_node] = nodeToDegree[u];
						nodeToIndex[end_node] = nodeToIndex[u];
					}
					bucket[nodeToDegree[u]].pop_back();
					nodeToDegree[u] = -1;
					nodeToIndex[u] = -1;
				}
				if (score[u] < max_w && score[u] > 0)
				{
					active[u] = true;
					bucket[score[u]].push_back(u);
					nodeToIndex[u] = bucket[score[u]].size() - 1;
					nodeToDegree[u] = score[u];
					nn = min(nn, nodeToDegree[u]);
				}
			}
		}
		if ((++cnt) * score[x] > 1000000)
		{
			cnt = 0;
		}
	}
	n_core = 0;
	m_core = 0;
	// for (int u = 0; u < n; ++u)
	// 	if (rank[u] == -1)
	// 	{
	// 		++n_core;
	// 		m_core += (int)E[u].size();
	// 	}
}

void CoreTree::load_graph(string bin_path)
{
	long long p = 0;
	FILE *fin = fopen((bin_path + "graph-dis.bin").c_str(), "rb");
	fread(&n, sizeof(int), 1, fin);
	fread(&m, sizeof(long long), 1, fin);
	deg = new int[n];
	dat = new int[m];
	con = new int *[n];
	nid = new int[n];
	old = new int[n];
	new_old = new int[n];
	fread(deg, sizeof(int), n, fin);
	fread(dat, sizeof(int), m, fin);
	fread(nid, sizeof(int), n, fin);
	fread(old, sizeof(int), n, fin);
	fclose(fin);
	for (int i = 0; i < n; ++i)
	{
		con[i] = dat + p;
		p += deg[i];
	}
}

bool CoreTree::get_edge(char *line, int &a, int &b, int num_cnt)
{
	if (!isdigit(line[0]))
		return false;
	vector<char *> v_num;
	int len = (int)strlen(line);
	for (int i = 0; i < len; ++i)
		if (!isdigit(line[i]) && line[i] != '.')
			line[i] = '\0';
		else if (i == 0 || !line[i - 1])
			v_num.push_back(line + i);
	if ((int)v_num.size() != num_cnt)
		return false;
	sscanf(v_num[0], "%d", &a);
	sscanf(v_num[1], "%d", &b);
	return true;
}

int CoreTree::get_num_cnt(string path)
{
	FILE *fin = fopen((path).c_str(), "r");
	char line[MAXLINE];
	int cnt = 0, min_cnt = 100;

	while (fgets(line, MAXLINE, fin) && cnt < 10)
	{
		if (!isdigit(line[0]))
			continue;
		vector<char *> v_num;
		int len = (int)strlen(line);
		for (int i = 0; i < len; ++i)
			if (!isdigit(line[i]) && line[i] != '.')
				line[i] = '\0';
			else if (i == 0 || !line[i - 1])
				v_num.push_back(line + i);
		if ((int)v_num.size() < 2)
			continue;
		min_cnt = min(min_cnt, (int)v_num.size());
		++cnt;
	}
	fclose(fin);
	return min_cnt;
}

void CoreTree::get_order(vector<int> *con, int n, int *o, int method)
{
	printf("method=%d\n", method);
	if (method == RANK_STATIC)
	{
		printf("Ranking Method = RANK_STATIC\n");
		// return val
		DV *f = new DV[n];
		for (int i = 0; i < n; ++i)
			f[i].id = i, f[i].val = con[i].size() * 1.0;
		sort(f, f + n);
		for (int i = 0; i < n; ++i)
			o[i] = f[i].id;
		delete[] f;
	}
}

void CoreTree::create_bin(string path, string graph_name, int rank_threads, int rank_method, int rank_max_minutes, int max_hops, bool merge_equv, string bin_path)
{
	string full_path = path + graph_name;
	cout << "full_path:" << full_path << endl;
	FILE *fin = fopen((full_path).c_str(), "r");
	char line[MAXLINE];
	int n = 0, a, b, num_cnt = get_num_cnt(full_path);
	vector<pair<int, int>> el;
	long long cnt = 0, m = 0;

	printf("Loading text, num_cnt=%d...\n", num_cnt);
	while (fgets(line, MAXLINE, fin))
	{
		if (!get_edge(line, a, b, num_cnt))
			continue;
		if (a < 0 || b < 0 || a == b)
			continue;
		el.push_back(make_pair(a, b));
		n = max(max(n, a + 1), b + 1);
		if ((++cnt) % (long long)10000000 == 0)
			printf("%lld lines finished\n", cnt);
	}
	fclose(fin);

	cout << "origional node:" << n << endl;
	vector<int> *con1 = new vector<int>[n];
	printf("Deduplicating...\n");

	for (size_t i = 0; i < el.size(); ++i)
	{
		con1[el[i].first].push_back(el[i].second);
		con1[el[i].second].push_back(el[i].first);
	}

	vector<int> degg;
	degg.resize(n);
	for (int i = 0; i < n; ++i)
		if (con1[i].size() > 0)
		{
			sort(con1[i].begin(), con1[i].end());
			int p = 1;
			for (int j = 1; j < (int)con1[i].size(); ++j)
				if (con1[i][j - 1] != con1[i][j])
					con1[i][p++] = con1[i][j];
			con1[i].resize(p);
			m += p;
			degg[i] = p;
		}
	// conn
	int id = 0;
	vector<Val> t;
	t.clear();
	for (int i = 0; i < n; ++i)
	{
		t.push_back(Val(i, degg[i]));
	}
	sort(t.begin(), t.end());
	id = t[0].id;

	queue<int> to_visit;
	vector<int> vistied;
	vistied.resize(n);
	vector<int> dis;
	dis.clear();
	dis.resize(n);

	vector<int> old_new;
	old_new.resize(n, -1);
	to_visit.push(id);
	vistied[id] = 1;
	cnt = 0;
	dis[id] = 0;
	old_new[id] = cnt++;
	while (!to_visit.empty())
	{
		int n_id = to_visit.front();
		to_visit.pop();

		for (vector<int>::iterator i = con1[n_id].begin();
			 i != con1[n_id].end();
			 i++)
		{
			if (vistied[*i] == 0)
			{
				vistied[*i] = 1;
				dis[*i] = dis[n_id] + 1;
				old_new[*i] = cnt++;
				to_visit.push(*i);
			}
		}
	}
	vector<pair<int, int>> n_el;
	n = 0;
	m = 0;

	for (int i = 0; i < el.size(); i++)
	{
		int a = el[i].first;
		int b = el[i].second;
		if (old_new[a] > -1 and old_new[b] > -1)
		{
			int n_a = old_new[a];
			int n_b = old_new[b];
			n_el.push_back(make_pair(n_a, n_b));
			n = max(max(n, n_a + 1), n_b + 1);
		}
	}

	printf("# conn nodes=%d\n", n);

	vector<int> *con = new vector<int>[n];

	for (long long i = 0; i < n_el.size(); ++i)
	{
		con[n_el[i].first].push_back(n_el[i].second);
		con[n_el[i].second].push_back(n_el[i].first);
	}
	for (int i = 0; i < n; ++i)
		if (con[i].size() > 0)
		{
			sort(con[i].begin(), con[i].end());
			int p = 1;
			for (int j = 1; j < (int)con[i].size(); ++j)
				if (con[i][j - 1] != con[i][j])
					con[i][p++] = con[i][j];
			con[i].resize(p);
			m += p;
			degg[i] = p;
		}
	printf("n=%d %d\n", n, m);
	long long *f1 = new long long[n];
	memset(f1, 0, sizeof(long long) * n);

	long long *f2 = new long long[n];
	memset(f2, 0, sizeof(long long) * n);

	double tt = omp_get_wtime();

	for (int i = 0; i < n; ++i)
		f1[i] = i, f2[i] = i;

	printf("Reordering...\n");
	int *f = new int[n];
	get_order(con, n, f, rank_method);

	int *old = new int[n], *nid = new int[n];
	for (int i = 0; i < n; ++i)
		nid[f[i]] = i;
	for (int i = 0; i < n; ++i)
		old[nid[i]] = i;

	for (int i = 0; i < n; ++i)
	{
		if (f1[i] != i)
			nid[i] = -nid[f1[i]] - 1;
		if (f2[i] != i)
			nid[i] = MAXN + nid[f2[i]];
	}
	cout << "n:" << n << ",m:" << m << endl;
	printf("Creating adjacency list...\n");
	int *dat = new int[m], *deg = new int[n], **adj = new int *[n];

	long long pos = 0;
	for (int i = 0; i < n; ++i)
	{
		adj[i] = dat + pos;
		pos += (int)con[old[i]].size();
	}
	memset(deg, 0, sizeof(int) * n);
	for (int i = 0; i < n; ++i)
	{
		int ii = old[i];
		for (int p = 0; p < (int)con[ii].size(); ++p)
		{
			int jj = con[ii][p];
			int j = nid[jj];
			adj[j][deg[j]++] = i;
		}
	}

	printf("Saving binary...\n");

	FILE *fout = fopen((bin_path + "graph-dis.bin").c_str(), "wb");
	fwrite(&n, sizeof(int), 1, fout);
	fwrite(&m, sizeof(long long), 1, fout);
	fwrite(deg, sizeof(int), n, fout);
	fwrite(dat, sizeof(int), m, fout);
	fwrite(nid, sizeof(int), n, fout);
	fwrite(old, sizeof(int), n, fout);
	fclose(fout);

	printf("Created binary file, n = %d, m = %lld\n", n, m);
}

void CoreTree::create_landmark()
{
	init(max_w);
	select_landmark();
	distance_landmark();
	update_bound();
	assign_group();
}
void CoreTree::select_landmark()
{
	landId.clear();

	vector<Val> l;
	for (int i = 0; i < n; ++i)
		l.push_back(Val(i, V[i].len));

	sort(l.begin(), l.end());

	int maxx = 99999, bs, id;
	vector<int> di;
	for (int i = 0; i < numLand; ++i)
		landId.push_back(l[i].id);
}

void CoreTree::distance_landmark()
{
	disLand.resize(numLand, vector<int>(n, 0));

	for (int i = 0; i < numLand; ++i)
	{
		int uid = landId[i];
		BFS(uid, disLand[i]);
	}
}

void CoreTree::update_bound()
{
	listLand.resize(numLand);
	for (int i = 0; i < numLand; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			listLand[i].push_back(Val(j, disLand[i][j]));
			update_lb(j, disLand[i][j]);
		}
		sort(listLand[i].begin(), listLand[i].end());
		int ecc = listLand[i][0].val;
		for (int j = 0; j < n; ++j)
		{
			update_lb(j, ecc - disLand[i][j]);
			update_ub(j, ecc + disLand[i][j]);
		}
	}
}

void CoreTree::assign_group()
{
	nodeGroup.resize(numLand);
	for (int j = 0; j < n; ++j)
	{
		int gid = 0, dis = POSINF;
		for (int k = 0; k < numLand; ++k)
		{
			if (disLand[k][j] < dis)
			{
				gid = k, dis = disLand[k][j];
			}
		}
		nodeGroup[gid].push_back(j);
	}
}

int CoreTree::get_unmch(int lid)
{
	int unmch = 0;
	for (int i = 0; i < nodeGroup[lid].size(); ++i)
	{
		int id = nodeGroup[lid][i];
		if (V[id].ub != V[id].lb)
		{

			unmch++;
		}
	}
	return unmch;
}

void CoreTree::update_lb(int v, int lb)
{
	if (lb > V[v].lb)
		V[v].lb = lb;
}

void CoreTree::update_ub(int v, int ub)
{
	if (ub < V[v].ub)
		V[v].ub = ub;
}

void CoreTree::compute_ecc()
{
	for (int lid = 0; lid < numLand; ++lid)
	{
		compute_ecc(lid);
	}
}

void CoreTree::txt_to_file(string graph_name)
{
	string file_name = "txt/" + graph_name;
	FILE *out = fopen(file_name.c_str(), "w");
	if (out == NULL)
	{
		fprintf(stderr, "can't write the graph file \n");
		cerr << graph_name << endl;
		exit(1);
	}

	vector<pair<int, int>> v;
	for (int i = 0; i < n; i++)
	{
		fprintf(out, "%d %d\n", i, V[i].lb);
	}

	fclose(out);
}

void CoreTree::compute_ecc(int lid)
{
	int lambcnt = -1, cnt = 0, lambda = -1;
	int core_num = 0, branch_num = 0;
	int max_id = -1, max_lb = -1, max_difference = 0;
	// int unmch = get_unmch(lid);
	vector<int> visited;
	visited.resize(n);
	vector<int> node_vis;
	vector<int> node_continue;
	node_vis.resize(n, 0);
	node_continue.resize(n, 0);
	record_dis.resize(n);
	used_dis.resize(n);

	while (true)
	{
		lambcnt++;
		int uid = listLand[lid][lambcnt].id;
		int group_index = node_g[uid];
		used_dis[group_mem[group_index].size()]++;
		lambda = listLand[lid][lambcnt].val;
		if (check(uid, visited))
		{
			Local_BFS(uid, lambda, lid);
			if (max_id != -1 && disLand[lid][max_id] + lambda <= max_lb)
			{
				break;
			}
			continue;
		}
		if (node_vis[uid])
		{
			continue;
		}
		vector<int> dis(n, 0);
		if (BFS_Judge(uid, node_vis, node_continue))
		{
			// cout << "uid:" << uid << ",deg:" << deg[uid] << endl;
			cnt++;

			int ecc = NEGINF;
			int group_index = node_g[uid];
			int left = group_l[group_index], right = group_r[group_index];
			max_difference = 0;
			if (deg[uid] == 2)
			{
				int ecc1 = NEGINF, ecc2 = NEGINF;
				for (int i = 0; i < n; i++)
				{
					ecc1 = max(ecc1, Dij_dis[left][i]);
					ecc2 = max(ecc2, Dij_dis[right][i]);
				}
				V[left].lb = V[left].ub = ecc1;
				V[right].lb = V[right].ub = ecc2;
				for (int i = 0; i < nodeGroup[lid].size(); ++i)
				{
					int id = nodeGroup[lid][i];
					int lb = V[id].lb, ub = V[id].ub;
					int u, left_distance, right_distance;
					if (node_g[id] == group_index || id == left || id == right)
					{
						left_distance = Dij_dis[left][id];
						right_distance = Dij_dis[right][id];
					}
					else
					{
						u = binarySearch(group_d[group_index], (group_len[group_index] + Dij_dis[right][id] - Dij_dis[left][id]) / 2, group_len[group_index], Dij_dis[left][id], Dij_dis[right][id]);
						left_distance = Dij_dis[left][id] + group_d[group_index][u];
						right_distance = Dij_dis[right][id] + group_len[group_index] - group_d[group_index][u];
					}
					if (lb < disLand[lid][id] + lambda)
					{
						int d = min(left_distance, right_distance);
						update_lb(id, d);
						if (V[id].lb == disLand[lid][id] + lambda)
						{
							continue;
						}
						if (disLand[lid][id] + lambda - V[id].lb > max_difference)
						{
							max_difference = disLand[lid][id] + lambda - V[id].lb;
							max_id = id;
							max_lb = V[id].lb;
						}
						d_l[left] = min(d_l[left], V[id].lb - Dij_dis[left][id]);
						d_l[right] = min(d_l[right], V[id].lb - Dij_dis[right][id]);
					}
				}
				// for (int i = 0;i < n;i++) {
				// 	int id = listLand[0][i].id;
				// 	cout << "id:" << id << ",lb:" << V[id].lb << ",ub:" << disLand[lid][id] + lambda << endl;
				// }

				// cout << "cnt:" << cnt << ",unmch:" << max_difference << endl;
				if (max_difference == 0)
					break;
			}
			else
			{
				for (int i = 0; i < n; i++)
				{
					ecc = max(ecc, Dij_dis[uid][i]);
				}
				V[uid].lb = V[uid].ub = ecc;
				for (int i = 0; i < nodeGroup[lid].size(); ++i)
				{
					int id = nodeGroup[lid][i];
					int lb = V[id].lb, ub = V[id].ub;
					if (lb < disLand[lid][id] + lambda)
					{
						update_lb(id, Dij_dis[uid][id]);
						if (V[id].lb == disLand[lid][id] + lambda)
						{
							continue;
						}
						if (disLand[lid][id] + lambda - V[id].lb > max_difference)
						{
							max_difference = disLand[lid][id] + lambda - V[id].lb;
							max_id = id;
							max_lb = V[id].lb;
						}
						d_l[uid] = min(d_l[uid], V[id].lb - Dij_dis[uid][id]);
					}
				}
				// for (int i = 0;i < n;i++) {
				// 	int id = listLand[0][i].id;
				// 	cout << "id:" << id << ",lb:" << V[id].lb << ",ub:" << disLand[lid][id] + lambda << endl;
				// }
				// cout << "cnt:" << cnt << ",unmch:" << max_difference << endl;
				if (max_difference == 0)
					break;
			}
		}
		else
		{
			// cout << "uid:" << uid << endl;
			if (max_id != -1 && disLand[lid][max_id] + lambda <= max_lb)
			{
				break;
			}
		}
	}
}

#endif /* CORETREELABELLING_H_ */
