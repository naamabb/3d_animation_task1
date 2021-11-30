#pragma once
#include "igl/opengl/glfw/Viewer.h"
#include "igl/aabb.h"

typedef std::set<std::pair<double, int> > PriorityQueue;

class SandBox : public igl::opengl::glfw::Viewer
{
public:

	SandBox();
	~SandBox();
	void Init(const std::string& config);
	// load list of meshes from txt file and
// and construct for each data structure for decimation.
	void ObjectLoader(const std::string config, bool is_shortest);
	//simplification function: input - number of faces to delete
	void SimplifyMesh(int n);
	void SetQueue(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	size_t getQsize(int id);

	

	//Q 10 - 11
	Eigen::Matrix4d* calculate_plane_err(Eigen::Vector3d normal, Eigen::Vector3d v);
	void update_vertex_errors(int v, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
 
	void calculate_cost_and_placement(const int e, double& cost, Eigen::RowVectorXd& p);
	Eigen::Matrix4d* derived_matrix(Eigen::Matrix4d& K);
	Eigen::Vector4d* min_choice(int v1, int v2, Eigen::Matrix4d& K);
	void print_data_structure();
	bool myCollapse_edge(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void kill_edge(int data_idx, int e);
	void _setQueue(Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	double doubleVariable;
	bool is_shortest;
private:
	// Prepare array-based edge data structures and priority queue
		// Prepare array-based edge data structures and priority queue
	std::vector<Eigen::VectorXi*> EMAP;//connects faces to edges
	//E.row(EMAP(f+i*F.rows)) --> [s,d], where E.row(i) --> [s,d]
	std::vector<Eigen::MatrixXi*> E, EF, EI;
	//E – edges <index of source vertex, index of destination vertex>
	//EF – connects edges to faces
	//EI – connects edge to vertex index in triangle (0,1,2)
	typedef std::set<std::pair<double, int> > PriorityQueue;
	std::vector<PriorityQueue*> Q;
	std::vector<std::vector<PriorityQueue::iterator>> Qit;
	//Q – priority queue containing cost for every edge
	// If an edge were collapsed, we'd collapse it to these points:
	std::vector<Eigen::MatrixXd*> C;//C – position of the new vertex after collapsing the corresponding edge (in model
	//coordinates).
	std::vector<int> num_collapsed;

	//Q 10 - 11
	
	std::vector<std::vector<Eigen::Matrix4d*>> vertex_errors;
	
	
	void Animate();
};

