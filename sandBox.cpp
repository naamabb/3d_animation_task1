#include "tutorial/sandBox/sandBox.h"
#include <igl/edge_flaps.h>
#include <igl/collapse_edge.h>
#include "Eigen/dense"
#include <functional>


#include <igl/circulation.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include <set>
#include <Eigen/LU>

#include "igl/circulation.h"




SandBox::SandBox()
{
	

}

void SandBox::Init(const std::string &config)
{
   
	std::string item_name;
	std::ifstream nameFileout;
	doubleVariable = 0;
	nameFileout.open(config);
	if (!nameFileout.is_open())
	{
		std::cout << "Can't open file "<<config << std::endl;
	}
	else
	{
		
		while (nameFileout >> item_name)
		{
			std::cout << "openning " << item_name << std::endl;
			load_mesh_from_file(item_name);
			
			parents.push_back(-1);
			data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
			data().show_overlay_depth = false;
			data().point_size = 10;
			data().line_width = 2;
			data().set_visible(false, 1);

			
		}
		nameFileout.close();
	}
	MyTranslate(Eigen::Vector3d(0, 0, -1), true);
	
	data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));

}

void SandBox::ObjectLoader(const std::string config, bool is_shortest)
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    std::string item_name;
    std::ifstream nameFileout;
    nameFileout.open(config);
    if (!nameFileout.is_open())
    {
        std::cout << "Can't open file " << config << std::endl;
    }
    else
    {

        while (nameFileout >> item_name)
        {
            // Create new data slot and set to selected
            std::cout << "openning " << item_name << std::endl;
            load_mesh_from_file(item_name);
            parents.push_back(-1);
            data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
            data().show_overlay_depth = false;
            data().point_size = 10;
            data().line_width = 2;
            data().set_visible(false, 1);
            num_collapsed.push_back(0);
    /**        if (!(data().F.rows() == 0 && data().V.rows() == 0))
            {
                append_mesh();
            }
            data().clear();
            MatrixXd V;
            MatrixXi F;
            read_triangle_mesh(item_name, V, F);
            data().set_mesh(V, F);*/
            if (is_shortest) {
                SetQueue(data().V, data().F);
            }
            else {
                _setQueue(data().V, data().F);

            }
            



        }
    }
}

void SandBox::SimplifyMesh(int n)
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    size_t curIdx = mesh_index(data().id);
    if (!Q[curIdx]->empty())
    {
        bool something_collapsed = false;
        // collapse edge
        const int max_iter = n;
        Eigen::MatrixXd nV = data().V;
        Eigen::MatrixXi nF = data().F;
        
        for (int j = 0; j < max_iter; j++)
        {

            if (is_shortest) {
                if (!collapse_edge(
                    shortest_edge_and_midpoint, nV, nF, *E[curIdx], *EMAP[curIdx], *EF[curIdx], *EI[curIdx], *Q[curIdx], Qit[curIdx], *C[curIdx]))
                {
                    break;
                }
            }
            else {
                if (!myCollapse_edge(nV, nF)) {
                    break;
                }
            }

            something_collapsed = true;
            num_collapsed[curIdx]++;
            printf("something collapsed: %d\n", num_collapsed[curIdx]);
        }

        if (something_collapsed)
        {
            data().clear();
            data().set_mesh(nV, nF);
            data().set_face_based(true);
            data().dirty = 157;
        }
    }
}


void SandBox::SetQueue(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    // Prepare array-based edge data structures and priority queue
    VectorXi* EMAPtmp = new VectorXi();
    MatrixXi* Etmp = new MatrixXi(), * EFtmp = new MatrixXi(), * EItmp = new MatrixXi();
    PriorityQueue *Qtmp = new PriorityQueue();
    std::vector<PriorityQueue::iterator> QitTmp;
    
    // If an edge were collapsed, we'd collapse it to these points:
    edge_flaps(F, *Etmp, *EMAPtmp, *EFtmp, *EItmp);
    C.push_back(new MatrixXd(Etmp->rows(), V.cols()));
    QitTmp.resize(Etmp->rows());
    VectorXd costs(Etmp->rows());
    Qtmp->clear();
    

    for (int e = 0; e < Etmp->rows(); e++)
    {
        double cost = e;
        RowVectorXd p(1, 3);
        shortest_edge_and_midpoint(e, V, F, *Etmp, *EMAPtmp, *EFtmp, *EItmp, cost, p);
        C.back()->row(e) = p;
        costs(e) = cost;
        QitTmp[e] = Qtmp->insert(std::pair<double, int>(cost, e)).first;
    }


    EMAP.push_back(EMAPtmp);
    E.push_back(Etmp);
    EF.push_back(EFtmp);
    EI.push_back(EItmp);
    Q.push_back(Qtmp);
    Qit.push_back(*new std::vector<PriorityQueue::iterator>(QitTmp));

    
}

void SandBox::_setQueue(Eigen::MatrixXd& V, Eigen::MatrixXi& F) 
{
    using namespace std;
    using namespace Eigen;
    using namespace igl;
    // Prepare array-based edge data structures and priority queue
    VectorXi* EMAPtmp = new VectorXi();
    MatrixXi* Etmp = new MatrixXi(), * EFtmp = new MatrixXi(), * EItmp = new MatrixXi();
    PriorityQueue* Qtmp = new PriorityQueue();
    std::vector<PriorityQueue::iterator> QitTmp;
    std::vector<Eigen::Matrix4d*>* errors = new std::vector<Eigen::Matrix4d*>;
    errors->resize(V.rows());
    vertex_errors.push_back(*errors);

    // If an edge were collapsed, we'd collapse it to these points:
    edge_flaps(F, *Etmp, *EMAPtmp, *EFtmp, *EItmp);
    C.push_back(new MatrixXd(Etmp->rows(), V.cols()));
    QitTmp.resize(Etmp->rows());
    VectorXd costs(Etmp->rows());
    Qtmp->clear();

    EMAP.push_back(EMAPtmp);
    E.push_back(Etmp);
    EF.push_back(EFtmp);
    EI.push_back(EItmp);

    for (int e = 0; e < Etmp->rows(); e++)
    {
        VectorXi edge = Etmp->row(e);
        update_vertex_errors(edge.x(),V,F);
        update_vertex_errors(edge.y(),V,F);

    }

    for (int e = 0; e < Etmp->rows(); e++)
    {
        double cost = e;
        RowVectorXd p(1, 3);
        calculate_cost_and_placement(e, cost, p);
        C.back()->row(e) = p;
        costs(e) = cost;
        QitTmp[e] = Qtmp->insert(std::pair<double, int>(cost, e)).first;
    }



    Q.push_back(Qtmp);
    Qit.push_back(*new std::vector<PriorityQueue::iterator>(QitTmp));




}


size_t SandBox::getQsize(int id)
{
    return Q[mesh_index(id)]->size();
}

Eigen::Matrix4d* SandBox::calculate_plane_err(Eigen::Vector3d normal, Eigen::Vector3d v)
{
    double d = v.transpose() * normal;
    d *= -1.0;
    Eigen::Vector4d p;
    p<< normal, d;
    Eigen::Matrix4d tmp = p * p.transpose();
    Eigen::Matrix4d* K = new Eigen::Matrix4d(tmp);
    
    return K;
}

void SandBox::update_vertex_errors(int v, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    int data_idx = mesh_index(data().id);
   
    Eigen::Matrix4d* K =  new Eigen::Matrix4d(Eigen::Matrix4d::Zero());
    Eigen::Matrix4d zero = Eigen::Matrix4d::Zero();
    Eigen::Vector3i face;

    for (int i = 0; i < F.rows(); i++)
    {

        face  =  F.row(i);
        if (face.x() == 0 && face.y() == 0 && face.z() == 0)
        {
            continue; // the face already collapsed
        }
        if (face.x() == v || face.y() == v || face.z() == v)
        {
            
            *K += *calculate_plane_err(data().F_normals.row(i).normalized(), data().V.row(v));

        }


    }
    vertex_errors[data_idx][v] = K;
}



Eigen::Matrix4d* SandBox::derived_matrix(Eigen::Matrix4d& K) 
{
    Eigen::Matrix4d* res = new Eigen::Matrix4d();
    res->row(0) = K.row(0);
    res->row(1) = K.row(1);
    res->row(1)[0] = K.row(0)[1];
    res->row(2) = K.col(2);
    res->row(2)[3] = K.row(2)[3];
    res->row(3) = *new Eigen::RowVector4d(0, 0, 0, 1);
    return res;

    
}

Eigen::Vector4d* SandBox::min_choice(int v1, int v2, Eigen::Matrix4d& K) 
{
    int idx = mesh_index(data().id);
    double cost1, cost2, cost3;
    Eigen::Vector3d V1 = data().V.row(v1), V2 = data().V.row(v2);
    Eigen::Vector4d t1,t2,t3,t;
    t1 << (V1 + V2) / 2, 1;
    t2 << V1, 1;
    t3 << V2, 1;
    if (t1.transpose() * K * t1 < t2.transpose() * K * t2)
    {
        t = t1;
    }
    else {
        t = t2;
    }
    if (t.transpose() * K * t > t3.transpose() * K * t3)
    {
        t = t3;
    }
    Eigen::Vector4d *res = new Eigen::Vector4d(t);

    return res;

}



void SandBox::calculate_cost_and_placement(const int e, double& cost, Eigen::RowVectorXd& p)
{
    using namespace Eigen;
    int curIdx = mesh_index(data().id);
    Vector2i edge = E[curIdx]->row(e);
    Matrix4d K = *vertex_errors[curIdx][edge.x()] + *vertex_errors[curIdx][edge.y()];
    Vector4d new_pos;   
    Matrix4d dK = *derived_matrix(K);
    if (dK.determinant() > 0) {
        new_pos = dK.inverse() * Vector4d(0, 0, 0, 1);
    }
    else {
        new_pos = *min_choice(edge.x(), edge.y(), K);
    }
        
   
    double *cost_res = new double (new_pos.transpose() * K * new_pos);
    if (new_pos[3] != 0) {
        new_pos[0] /= new_pos[3];
        new_pos[1] /= new_pos[3];
        new_pos[2] /= new_pos[3];
    }
    p = *new Vector3d(new_pos.head(3));
    cost = *cost_res;
}

void SandBox::print_data_structure()
{
    size_t curIdx = mesh_index(data().id);
    std::cout << "V:\n" << data().V << std::endl;
    std::cout << "F:\n" << data().F << std::endl;
    std::cout << "E:\n" << *E[curIdx] << std::endl;
    std::cout << "EF:\n" << *EF[curIdx] << std::endl;
    std::cout << "EI:\n" << *EI[curIdx] << std::endl;
    std::cout << "EMAP:\n" << *EMAP[curIdx] << std::endl;
    std::cout << "F_normals:\n" << data().F_normals << std::endl;
}

#define COLLAPSE_EDGE_NULL 0;

void SandBox::kill_edge(int data_idx ,int e)
{
    (*E[data_idx])(e,0) = COLLAPSE_EDGE_NULL;
    (*E[data_idx])(e,1) = COLLAPSE_EDGE_NULL;
    (*EF[data_idx])(e,0) = COLLAPSE_EDGE_NULL;
    (*EF[data_idx])(e,1) = COLLAPSE_EDGE_NULL;
    (*EI[data_idx])(e,0) = COLLAPSE_EDGE_NULL;
    (*EI[data_idx])(e,1) = COLLAPSE_EDGE_NULL;

}

bool SandBox::myCollapse_edge(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    int curIdx = mesh_index(data().id);
    std::vector<Eigen::Matrix4d*> v_errors = vertex_errors[curIdx];

    if (Q[curIdx]->empty()) {
        return false;
    }

    std::pair<double, int> p = *(*Q[curIdx]).begin();
    int e = p.second;
    double cost = p.first;
    Eigen::Vector3d* n_pos = new Eigen::Vector3d(C[curIdx]->row(e));
   // std::cout << "edge to remove: " << e << std::endl;


        const auto& derived_matrix = [](Eigen::Matrix4d& K, Eigen::Matrix4d& res)
        {
       // Eigen::Matrix4d* res = new Eigen::Matrix4d();
        res.row(0) = K.row(0);
        res.row(1) = K.row(1);
        res.row(1)[0] = K.row(0)[1];
        res.row(2) = K.col(2);
        res.row(2)[3] = K.row(2)[3];
        res.row(3) = *new Eigen::RowVector4d(0, 0, 0, 1);
        

        };
        const auto& min_choice = [&V](int v1, int v2, Eigen::Matrix4d& K, Eigen::Vector4d& res)
        {
            Eigen::Vector3d V1 = V.row(v1), V2 = V.row(v2);
            Eigen::Vector4d t1, t2, t3, t;
            t1 << (V1 + V2) / 2, 1;
            t2 << V1, 1;
            t3 << V2, 1;
            if (t1.transpose() * K * t1 < t2.transpose() * K * t2)
            {
                t = t1;
            }
            else {
                t = t2;
            }
            if (t.transpose() * K * t > t3.transpose() * K * t3)
            {
                t = t3;
            }
            res = *new Eigen::Vector4d(t);

        

        };
        const auto& calculate_cost_and_placement = [&derived_matrix, &min_choice, &v_errors](
            const int e,
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& /*F*/,
            const Eigen::MatrixXi& E,
            const Eigen::VectorXi& /*EMAP*/,
            const Eigen::MatrixXi& /*EF*/,
            const Eigen::MatrixXi& /*EI*/,
            double& cost,
            Eigen::RowVectorXd& p)
        {
            using namespace Eigen;
            Vector2i edge = E.row(e);
            Matrix4d K = *v_errors[edge.x()] + *v_errors[edge.y()];
            Vector4d new_pos;
            Matrix4d dK;
            derived_matrix(K,dK);
            if (dK.determinant() > 0) {
                std::cout << "inverse true" << std::endl;
                new_pos = dK.inverse() * Vector4d(0, 0, 0, 1);
            }
            else {
                min_choice(edge.x(), edge.y(), K, new_pos);
            }


            double* cost_res = new double(new_pos.transpose() * K * new_pos);
            if (new_pos[3] != 0) {
                new_pos[0] /= new_pos[3];
                new_pos[1] /= new_pos[3];
                new_pos[2] /= new_pos[3];
            }
            p = *new Vector3d(new_pos.head(3));
            cost = *cost_res;

        };
        bool res = igl::collapse_edge(calculate_cost_and_placement, V, F, *E[curIdx], *EMAP[curIdx], *EF[curIdx], *EI[curIdx], *Q[curIdx], Qit[curIdx], *C[curIdx]);
        if (res) {
            std::cout << "edge " << e << " cost = " << cost << " new v position " << n_pos->x() << " " << n_pos->y() << " " << n_pos->z() << std::endl;
        }
        return res;
        

}

/*bool SandBox::myCollapse_edge(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    int data_idx = mesh_index(data().id);
    if (Q[data_idx]->empty()) {
        return false;
    }

    std::pair<double, int> p = *(*Q[data_idx]).begin();
    int e = p.second;
    double cost = p.first;
    Q[data_idx]->erase(Q[data_idx]->begin());
    std::cout << "edge to remove: " << e << std::endl;

    Eigen::Vector3d* n_pos = new Eigen::Vector3d(C[data_idx]->row(e));
    const int eflip = (*E[data_idx])(e,0) > (*E[data_idx])(e,1);
    // source and destination
    const int s = eflip ? (*E[data_idx])(e,1) : (*E[data_idx])(e,0);
    const int d = eflip ? (*E[data_idx])(e,0) : (*E[data_idx])(e,1);
    //grab neighbors of d
    const std::vector<int> nV2Fd = igl::circulation(e, !eflip, *EMAP[data_idx], *EF[data_idx], *EI[data_idx]);
    // move source and destination to new position
   
    V.row(s) = *n_pos;
    V.row(d) = *n_pos;

    // update edge info
 // for each flap
    const int m = F.rows();
    for (int side = 0; side < 2; side++)
    {
        const int f = (*EF[data_idx])(e,side);
        const int v = (*EI[data_idx])(e,side);
        const int sign = (eflip == 0 ? 1 : -1) * (1 - 2 * side);
        // next edge emanating from d
        const int e1 = (*EMAP[data_idx])(f + m * ((v + sign * 1 + 3) % 3));
        // prev edge pointing to s
        const int e2 = (*EMAP[data_idx])(f + m * ((v + sign * 2 + 3) % 3));
        // face adjacent to f on e1, also incident on d
        const bool flip1 = (*EF[data_idx])(e1,1) == f;
        const int f1 = flip1 ? (*EF[data_idx])(e1,0) : (*EF[data_idx])(e1,1);

        // across from which vertex of f1 does e1 appear?
        const int v1 = flip1 ? (*EI[data_idx])(e1,0) : (*EI[data_idx])(e1,1);
        // Kill e1
        kill_edge(data_idx, e1);
        // Kill f
        F(f, 0) = COLLAPSE_EDGE_NULL;
        F(f, 1) = COLLAPSE_EDGE_NULL;
        F(f, 2) = COLLAPSE_EDGE_NULL;
        // map f1's edge on e1 to e2

        (*EMAP[data_idx])(f1 + m * v1) = e2;
        // side opposite f2, the face adjacent to f on e2, also incident on s
        const int opp2 = ((*EF[data_idx])(e2,0) == f ? 0 : 1);

        (*EF[data_idx])(e2,opp2) = f1;
        (*EI[data_idx])(e2,opp2) = v1;
        // remap e2 from d to s
        (*E[data_idx])(e2,0) = (*E[data_idx])(e2,0) == d ? s : (*E[data_idx])(e2,0);
        (*E[data_idx])(e2,1) = (*E[data_idx])(e2,1) == d ? s : (*E[data_idx])(e2,1);
        update_vertex_errors(s,V,F);
        //update_vertex_errors(d,V,F);
        for (auto f : nV2Fd)
        {
            for (int v = 0; v < 3; v++)
            {
                if (F(f, v) == d)
                {
                    const int flip1 = ((*EF[data_idx])((*EMAP[data_idx])(f + m * ((v + 1) % 3)), 0) == f) ? 1 : 0;
                    const int flip2 = ((*EF[data_idx])((*EMAP[data_idx])(f + m * ((v + 2) % 3)), 0) == f) ? 0 : 1;

                    (*E[data_idx])((*EMAP[data_idx])(f + m * ((v + 1) % 3)), flip1) = s;
                    (*E[data_idx])((*EMAP[data_idx])(f + m * ((v + 2) % 3)), flip2) = s;
                    F(f, v) = s;
                    break;
                }
            }
            for (int i = 0; i < 3; i++)
            {
                int cur_edge = (*EMAP[data_idx])(f + i * F.rows());
                if (cur_edge != e) {
                    double _cost = cur_edge;
                    Eigen::RowVectorXd _p(1, 3);
                    calculate_cost_and_placement(cur_edge, _cost, _p);
                    C[data_idx]->row(cur_edge) = _p;
                    Q[data_idx]->erase(Qit[data_idx][cur_edge]);
                    Qit[data_idx][cur_edge] = Q[data_idx]->insert(std::pair<double, int>(_cost, cur_edge)).first;

                }


            }
        }
        // Finally, "remove" this edge and its information
        kill_edge(data_idx, e);
        std::cout << "edge " << e << " cost = " << cost << " new v position " << n_pos->x() << " " << n_pos->y()<<" " << n_pos->z() << std::endl;
        return true;






    }
}*/

    


SandBox::~SandBox()
{

}

void SandBox::Animate()
{
	if (isActive)
	{
		
		
		
	}
}


