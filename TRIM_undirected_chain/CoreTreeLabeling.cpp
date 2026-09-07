#include "CoreTreeLabelling.h"
// #include "IFECC.h"

#include <map>

void query_dis(string path, string graph_name, int max_w, int numLand, string bin_path)
{
	CoreTree ct(path);
	ct.load_label(max_w, bin_path, numLand);
	double tm = get_time();
	ct.create_landmark();
	ct.Compute_group();
	ct.Merge_root_dis();
	ct.compute_ecc();
	int max_dis = 0, min_dis = INT_MAX;
	int r = 10000;
	int d = 0;
	for (int i = 0; i < ct.n; ++i)
	{
		r = min(r, ct.V[i].lb);
		d = max(d, ct.V[i].lb);
	}
	cout << "r = " << r << ", d = " << d << endl;
	printf("%f\n", get_time() - tm);
	ct.txt_to_file(graph_name);
}

int main(int argc, char *argv[])
{
	setvbuf(stdout, NULL, _IONBF, 0);
	setvbuf(stderr, NULL, _IONBF, 0);
	string bin_path = "fb/";
	if (argc > 1)
	{
		if (strcmp(argv[1], "txt-to-bin") == 0)
		{
			CoreTree::create_bin(argv[2], argv[3],						 /*path*/
								 argc > 4 ? atoi(argv[4]) : 1,			 /*rank_threads*/
								 argc > 5 ? atoi(argv[5]) : RANK_STATIC, /*rank_method*/
								 argc > 6 ? atoi(argv[6]) : 60,			 /*rank_max_minutes*/
								 argc > 7 ? atoi(argv[7]) : 3,			 /*max_hops*/
								 argc > 8 ? atoi(argv[7]) : false,
								 bin_path);
		}
		else if (strcmp(argv[1], "decompose_tree") == 0)
		{
			CoreTree ct(argv[2]);							/*path*/
			ct.decompose_tree(atoi(argv[3]),				/*max_w*/
							  argc > 4 ? atoi(argv[4]) : 1, /*n_threads*/
							  bin_path);
		}
		else if (strcmp(argv[1], "decompose_bt") == 0)
		{
			CoreTree ct(argv[2]);							/*path*/
			ct.decompose_tree(atoi(argv[3]),				/*max_w*/
							  argc > 4 ? atoi(argv[4]) : 1, /*n_threads*/
							  bin_path);
			ct.save_tmp_graph(atoi(argv[3]), bin_path);
		}
		else if (strcmp(argv[1], "query-dis") == 0)
			query_dis(argv[2], /*lambmark*/
					  argv[3],
					  atoi(argv[4]), /*max_w*/
					  atoi(argv[5]), /*n_pairs*/
					  bin_path);
	}

	return 0;
}
