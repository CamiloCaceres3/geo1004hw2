#include <iostream>
#include <list>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stack>
#include <math.h>

#include <unordered_map>

#include "DCEL.hpp"

// forward declarations; these functions are given below main()
void DemoDCEL();
void printDCEL(DCEL & D);


/* 
  Example functions that you could implement. But you are 
  free to organise/modify the code however you want.
  After each function you should have a DCEL without invalid elements!
*/
// 1.
void importOBJ(DCEL & D, const char *file_in ,   std::unordered_map< HalfEdge*, std::vector<int> > &hemap, std::unordered_map<int, Vertex*> &umap) {
  std::string  input =  file_in;
  input = "../../" + input;
  std::cout << "Reading file: " << input << std::endl;
  std::ifstream infile(input.c_str(), std::ifstream::in);
  if (!infile)
  {
    std::cerr << "Input file not found.\n";
  }
  std::string cursor;
  std::string line = "";
  int id_vertex =1;
  int afaces = 0;
  int a = 0;
  while (  std::getline(infile, line))
  {
    std::istringstream linestream(line);
    linestream >> cursor;
    double x,y,z;
    int h,j,k;
    if(cursor == "v")
    {
      for(int i = 0; i < 3; i++)
      {
        linestream >> cursor;
        if(i == 0)
        {
          x = std::stod(cursor);
        }
        else if(i==1)
        {
          y = std::stod(cursor);
        }
        else if( i ==2)
        {
          z = std:: stod(cursor);
        }
      }
      // for each v one vertex
      Vertex* v0 = D.createVertex(x,y,z);
      umap.insert({id_vertex, v0});
      id_vertex  ++;
    }
    else if(cursor =="f")
    {
      for(int i = 0; i < 3; i++)
      {
        linestream >> cursor;
        if(i == 0)
        {
          h = std::stoi(cursor);
        }
        else if(i==1)
        {
          j = std::stoi(cursor);
        }
        else if( i ==2)
        {
          k = std::stoi(cursor);
        }
      }
      HalfEdge* e0 = D.createHalfEdge();
      HalfEdge* e1 = D.createHalfEdge();
      HalfEdge* e2 = D.createHalfEdge();
      
      Face* f0 = D.createFace();
      f0->exteriorEdge = e0;

      e0->origin = umap[h];
      e0->destination = umap[j];
      e0->next = e1;
      e0->prev = e2;
      e0->incidentFace = f0;
      std::vector<int> v {h,j};
      hemap.insert({e0,v});

      e1->origin = umap[j];
      e1->destination = umap[k];
      e1->next = e2;
      e1->prev = e0;
      e1->incidentFace = f0;
      std::vector<int> v1 {j,k};
      hemap.insert({ e1, v1});

      e2->origin = umap[k];
      e2->destination = umap[h];
      e2->next = e0;
      e2->prev = e1;
      e2->incidentFace = f0;
      std::vector<int> v2 {k,h};
      hemap.insert({ e2, v2});
    } 
  a++;
  //for each f 6 edges and for now 1 face (interior)
  }
  for(auto ed: hemap)
  {
    for(auto ed1: hemap)
      if((ed.second[0] == ed1.second[1] && ed.second[1] == ed1.second[0]&& ed.first != ed1.first) || (ed.second[0] == ed1.second[0] && ed.second[1] == ed1.second[1] && ed.first != ed1.first)  )
      {
        ed.first->twin = ed1.first;
        ed1.first-> twin = ed.first;
        //hemap.erase(ed.first);
        //hemap.erase(ed1.first);
      }
  }
  // auto it = hemap.begin()->first;
  //D.infiniteFace()->holes.push_back(it);
  std::cout << "hemap.size() is " << hemap.size() << "umap.size() is " << umap.size() <<std::endl;


  //printDCEL(D);

}
// 2.
void groupTriangles(DCEL & D, std::unordered_map< HalfEdge*, std::vector<int> > &hemap,  std::unordered_map< Face*, int> &facemap) {
  // to do
  //std::vector<int> meshes;
  // create a hashmap with the halfedge and an integer
  std::unordered_map< HalfEdge*, int> meshmap;
  for( auto & e : hemap)
  {
    meshmap.insert({e.first,0});
  }
  // buildings in .obj start at 1
  int build = 1;
  bool list_em = false;


  while(list_em == false)
  {
    //get a halfedge that has not been assigned to a mesh
    HalfEdge* c;
    std::unordered_map<HalfEdge*, int>::iterator it = meshmap.begin();
    while (it != meshmap.end()) {
        if (it->second == 0) {
            c = it->first;
            it == meshmap.end();
        } 
        it++;
    }
    //start the stack to go through all edges of a mesh
    std::stack<HalfEdge*> s;
    s.push(c);
    // store the edge into the list of holes of the infinite face
    D.infiniteFace()->holes.push_back(c);
    //traverse all the faces to get all meshes
    while(!s.empty())
    {
      HalfEdge* e = s.top();
      //std::cout << &e << std::endl;
      s.pop();
      facemap[e->incidentFace] = build;
      HalfEdge* e_start = e;
        do {
            meshmap[e] = 1;
            if(meshmap[e->twin] == 0)
              {
                s.push(e->twin);
              }
            e = e->next;
          } while ( e_start!=e) ; 
    }
    build ++;
    int ir =0;
    // you stop when you traversed all the edges
    for(auto c: meshmap)
    {
      if(c.second ==0)
      {
        ir ++;
      }
    }
    // stop the loop, you have found all meshes
    if(ir == 0)
    {
      list_em = true;
    }
  }
}

//std::cout << D.infiniteFace()->holes.size() << std::endl;

std::vector<double> cross(std::vector<double> V1, std::vector<double> V2)
{

  std::vector<double> v;
  v.push_back(V1[1] * V2[2] - V1[2] * V2[1]);
  v.push_back(-(V1[0] * V2[2] - V1[2] * V2[0]));
  v.push_back(V1[0] * V2[1] - V1[1] * V2[0]);
  return v;  
}

double dot(std::vector<double> V1, std::vector<double> V2)
{
  double d = 0.0;
  for (int i = 0; i < V1.size(); i++)
  {
    d = d + (V1[i] * V2[i]);
  }
  return d;
}


std::vector<double> getNormalVec(Vertex* &v1, Vertex* &v2, Vertex* &v3){
  //(v2-v1)X(v3-v1)
  std::vector<double> v12 {v1->x - v2->x, v1->y - v2->y, v1->z - v2->z};
  std::vector<double> v13 {v1->x - v3->x, v1->y - v3->y, v1->z - v3->z};
  std::vector<double> res = cross(v12,v13);
  return res;
}

std::vector<double> getPlaneCenter(Vertex* &v1, Vertex* &v2, Vertex* &v3){
  std::vector<double> vc {(v1->x + v2->x + v3->x)/3, (v1->y + v2->y + v3->y)/3, (v1->z + v2->z + v3->z)/3};
  return vc;
}

std::vector<double> vertToVec(Vertex* &v1){
  std::vector<double> vec1{v1->x,v1->y,v1->z};
  return vec1;
}

double signedVolume(std::vector<double> v1, std::vector<double> v2, std::vector<double> v3, std::vector<double> v4){
  std::vector<double> v1v4 {v1[0] - v4[0], v1[1] - v4[1], v1[2] - v4[2]};
  std::vector<double> v2v4 {v2[0] - v4[0], v2[1] - v4[1], v2[2] - v4[2]};
  std::vector<double> v3v4 {v3[0] - v4[0], v3[1] - v4[1], v3[2] - v4[2]};
  std::vector<double> v2v4v3v4 = cross(v2v4, v3v4);
  double dotv1v2v3v4 = dot(v1v4, v2v4v3v4);
  double res = dotv1v2v3v4/6;
  return res;
}


bool intersectCheck(Vertex* &v1, Vertex* &v2, Vertex* &v3, std::vector<double> origin, std::vector<double> destination){
  std::vector<double> vec1 = vertToVec(v1);
  std::vector<double> vec2 = vertToVec(v2);
  std::vector<double> vec3 = vertToVec(v3);
  //std::cout<<"origin "<< origin[0]<<" "<<  origin[1]<<" "<< origin[2] <<std::endl;
  //std::cout<<"destination "<< destination[0] <<" "<<destination[1]<<" "<< destination[2] <<std::endl;

  //1
  double vOrigin = signedVolume(vec1, vec2, vec3, origin);
  double vDest = signedVolume(vec1, vec2, vec3, destination);
  bool check1 = false;
  if ( vOrigin * vDest <= 0){
    check1 = true;
  }
  //2
  double v12 = signedVolume(origin, destination, vec1, vec2);
  double v13 = signedVolume(origin, destination, vec1, vec3);
  bool check2 = false;
  if ( v12 * v13 <= 0){
    check2 = true;
  }
//3
  double v21 = signedVolume(origin, destination, vec2, vec1);
  double v23 = signedVolume(origin, destination, vec2, vec3);
  bool check3 = false;
  if ( v21 * v23 <= 0){
    check3 = true;
  }
  //4
  double v31 = signedVolume(origin, destination, vec3, vec1);
  double v32 = signedVolume(origin, destination, vec3, vec2);
  bool check4 = false;
  if ( v31 * v32 <= 0){
    check4 = true;
  }
  if (check1 && check2 && check3 && check4){
    //std::cout<< v1->x<< " "<< v1->y<<" "<< v1->z<< std::endl;
    //std::cout<< v2->x<< " "<< v2->y<<" "<< v2->z<< std::endl;
    //std::cout<< v3->x<< " "<< v3->y<<" "<< v3->z<< std::endl;

    return true;

  }
  else{
    return false;
  };
}
void changeOrientation(DCEL & D, std::unordered_map< Face*, int> facemap, Face* f){
  Vertex* v_0 = f->exteriorEdge->origin;
  Vertex* v_1 = f->exteriorEdge->destination;
  Vertex* v_2 = f->exteriorEdge->next->destination;
  HalfEdge* e_0 = f->exteriorEdge;
  HalfEdge* e_1 = f->exteriorEdge->next;
  HalfEdge* e_2 = f->exteriorEdge->prev;

  e_0->origin = v_1;
  e_0->destination = v_0;
  e_0->next = e_2;
  e_0->prev = e_1;
  e_1->origin = v_2;
  e_1->destination = v_1;
  e_1->next = e_0;
  e_1->prev = e_2;
  e_2->origin = v_0;
  e_2->destination = v_2;
  e_2->next = e_1;
  e_2->prev = e_0;
  // take all edges from a triangle face and reverse them
}


void regGrowingOrientation(DCEL & D, std::unordered_map<Face*, int> fmap, std::unordered_map< Face*, int> facemap, HalfEdge*&ee , int count){
  //fmap[e->incidentFace] == 0;

  std::stack<HalfEdge*> s;
  s.push(ee);
  while(!s.empty() )
    {
      //take the first element of the stack
      HalfEdge* e = s.top();
      s.pop();

        if (e->origin == e->twin->origin ){
        changeOrientation(D, facemap, e->twin->incidentFace);
        }
        if (e->next->origin == e->next->twin->origin ){
        changeOrientation(D, facemap, e->next->twin->incidentFace);
        }
        if (e->prev->origin == e->prev->twin->origin ){
        changeOrientation(D, facemap, e->prev->twin->incidentFace);
        }
      
      fmap[e->incidentFace] = 1;
      if (fmap[e->next->twin->incidentFace] == 0){
          s.push(e->next->twin);
        };
      if (fmap[e->prev->twin->incidentFace] == 0){
          s.push(e->prev->twin);
          }
        if (fmap[e->twin->incidentFace] == 0){
          s.push(e->twin);
          }
        
        
  }

  //original edge, check twin orientation, 
  //                if orientation bad, change orientation, and perform reggrowing on twin->next
  //            check next edge twin orientation
  //                   if orientation bad, change orientation, and perform reggrowing on twin->next
  //            check last edge twin orientation
  //                if oriantation bad, change orientation, perform reggrowing on twin->next
  //recurse
  }

// 3.
void orientFace(DCEL & D, std::unordered_map< Face*, int> fmap, std::unordered_map< Face*, int> facemap, HalfEdge* e) {
    // to do
    //create normal & attach to center plane (normal vec * 10000 + point of face)
  std::cout<<"orientface"<<std::endl;
    std::vector<double> normaldir = getNormalVec(e->origin, e->destination, e->next->destination);
    std::vector<double> originnorm = getPlaneCenter(e->origin, e-> destination, e->next->destination);
    double scaledir = 10000;
    //std::vector<double> destnorm  { normaldir[0]*scaledir + originnorm[0], normaldir[1]*scaledir + originnorm[1], normaldir[2]*scaledir + originnorm[2]};
    std::vector<double> destnorm;
    destnorm.push_back(normaldir[0]*scaledir + originnorm[0]);
    destnorm.push_back(normaldir[1]*scaledir + originnorm[1]);
    destnorm.push_back(normaldir[2]*scaledir + originnorm[2]);

    float res = 0;
    int count = 0;
    for (auto & f2 : D.faces()){
      if ( e->incidentFace != f2->exteriorEdge->incidentFace){
        if (intersectCheck(f2->exteriorEdge->origin, f2->exteriorEdge->destination, f2->exteriorEdge->next->destination, originnorm, destnorm)){
          if (intersectCheck(f2->exteriorEdge->origin, f2->exteriorEdge->destination, f2->exteriorEdge->destination, originnorm, destnorm) ||
                  intersectCheck(f2->exteriorEdge->origin, f2->exteriorEdge->next->destination, f2->exteriorEdge->next->destination, originnorm, destnorm) ||
                  intersectCheck(f2->exteriorEdge->next->destination, f2->exteriorEdge->destination, f2->exteriorEdge->destination, originnorm, destnorm)){
            res = res + 0.5;
            std::cout<<"res "<<res<<std::endl;
          }
          count = count + 1;
          std::cout<<"count "<<count<<std::endl;
        }
      }
    }
    count = count - res;
              std::cout<<"count "<< count<<std::endl;

      if ( count % 2 != 0 && count != 0){
          std::cout<<"not even"<<std::endl;
          changeOrientation(D, facemap, e->incidentFace);
          orientFace(D, fmap, facemap, e);
          regGrowingOrientation(D, fmap, facemap, e, count);
      }
      else{
         std::cout<<"even"<<std::endl;
         regGrowingOrientation(D, fmap, facemap, e, count);
      }
      
    }



void OrientationOfMap(DCEL & D, std::unordered_map< Face*, int> fmap, std::unordered_map< Face*, int> facemap){
  for(  auto hedge: D.infiniteFace()->holes)
  {
    orientFace(D, fmap, facemap, hedge);
  }
}
// 4.




bool coplanar(HalfEdge* e1, HalfEdge* e2)
{

  std::vector<double> V1 {e1->destination->x - e1->next->destination->x, e1->destination->y - e1->next->destination->y , e1->destination->z - e1->next->destination->z };
  std::vector<double> V2 {e1->prev->origin->x - e1->origin->x, e1->prev->origin->y - e1->origin->y , e1->prev->origin->z - e1->origin->z };

  std::vector<double> V3 {e2->destination->x - e2->next->destination->x, e2->destination->y - e2->next->destination->y , e2->destination->z - e2->next->destination->z };
  std::vector<double> V4 {e2->prev->origin->x - e2->origin->x, e2->prev->origin->y - e2->origin->y , e2->prev->origin->z - e2->origin->z };



  std::vector<double> v  = cross(V1, V2);
  std::vector<double> v1  = cross(V3, V4);
  double d =  dot(v, V3);
  double d1 =  dot(v1, V2);
  double d2 =  dot(v, V4);
  double d3 =  dot(v1, V1);

  double D = -v[0]*e1->prev->origin->x  -v[1]*e1->prev->origin->y -v[2]*e1->prev->origin->z;

  double m1= sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
  double m2= sqrt(pow(v1[0],2)+pow(v1[1],2)+pow(v1[2],2));

  double angle1 = acos(abs(dot(v1,v)) /(m1*m2));

  double w = v[0]*e2->prev->origin->x  +v[1]*e2->prev->origin->y +v[2]*e2->prev->origin->z;

  double dist = abs(w+D) / m1; 

  if(  d < 0.01 && d1< 0.01 && d2 <0.01 && d3< 0.01)
  {
    return true;
  }
  else return false;
}

bool detect_hole(HalfEdge* e1, HalfEdge* e2)
{
  double xmax = e1->origin->x,ymax = e1->origin->y, xmin = e1->origin->x, ymin = e1->origin->y;
  //get bounding box
  const HalfEdge* e_start = e1;
        do {
            if(e1->origin->x > xmax) xmax = e1->origin->x;
            if(e1->origin->y > ymax) ymax = e1->origin->y;
            if(e1->origin->x < xmin) xmin = e1->origin->x;
            if(e1->origin->y < ymin) ymin = e1->origin->y;
            e1 = e1->next;
          } while ( e_start!=e1) ; 
            if(e2->origin->x < xmax && e2->origin->y < ymax && e2->origin->x > xmin &&  e2->origin->y > ymin ) return true; 
  return false  ;
}

bool is_hole(HalfEdge* e)
{
  Face* f = e->incidentFace;
  if(f->holes.empty())
  {
    return false;
  }
  else
  {
    for(auto  holes: f->holes)
    {
        const HalfEdge* e_start = holes;
        do {
            if(holes == e) return true;
            holes = holes->next;
            } while ( e_start!=holes) ; 
    }    
  }
  return false;

}


void mergeCoPlanarFaces(DCEL & D,   std::unordered_map< Face*, int> &facemap, std::unordered_map< HalfEdge*, int> &groupmap) {
  // to do
  //for each mesh of the model we merge faces
  for(  auto hedge: D.infiniteFace()->holes)
  {
    if(facemap[hedge->incidentFace] >=0)
    {
    std::cout << "Merging mesh:" <<  facemap[hedge->incidentFace]  << std::endl;
    //start with the mesh assigned in the list of holes of the infinite face
    HalfEdge* erg = hedge;
    // start a stack of  edges were each will be looked for coplanar faces
    std::stack<HalfEdge*> s;
    s.push(erg);
    while(!s.empty() )
    {
      //take the first element of the stack
      HalfEdge* e = s.top();
      s.pop();
      //if the edge has not been traversed
      if(groupmap[e]== 0)
      {
      //take 4 points and check if they are coplanar
      HalfEdge* te = e->twin;
      HalfEdge* n = e->next;
      HalfEdge* tn = te->next;
      HalfEdge* p = e->prev;    
      HalfEdge* tp = te->prev;
      Vertex* v1 = e->origin;
      Vertex* v2 = e->destination;
      Vertex* v3 = n->destination;
      Vertex* v4 = tn -> destination;
      Face* f = e->incidentFace;
      Face* tf = tn->incidentFace;

      //check if the are coplanar (all points should be coplanar)
      bool cop = true; 
      HalfEdge* tnn  = tn;
      HalfEdge* e_start = tnn;
      do {
          HalfEdge* ef  = e;
          HalfEdge* e_start1 = ef;
          do {
              if(!coplanar(tnn,ef)) 
              { cop = false;
                ef= e_start1;
              }
              else ef = ef->next;
              } while ( e_start1!=ef) ; 
          if(cop ==false)  tnn= e_start;
          else tnn = tnn->next;
        } while ( e_start!=tnn) ;  
      if(cop == true && (tp != n->twin &&   p != tn->twin ))
      {
        n->prev = tp;
        p->next = tn;
        tn->prev = p;
        tp->next = n;
        //delete the edge and the twin
        e -> eliminate();
        te -> eliminate();
        int hole_c = 0;
        //if the faces of the twin is the same -> we have a hole!
        //std::cout << "coplanar" <<  s.size()  << std::endl;
        if(tf == f)
        {
          // to check which edges belong to the hole
          //std::cout << "hole" <<  s.size()  << std::endl;
          if(detect_hole(n,f->exteriorEdge))
            {
              hole_c = 1;
             // std::cout << "n mayor p" <<  s.size()  << std::endl;
              f->holes.push_back(tn);           
            }
          else 
            {
              hole_c = 2;
              //std::cout << "p mayor n" <<  s.size()  << std::endl;
              f->holes.push_back(n);
            }
        }
        // if the faces are not the same -> delete the face
        else
        {
          for ( auto &g: tf-> holes)
          { 
            f->holes.push_back(g);
            const HalfEdge* e_start = g;
              do {
            g->incidentFace = f;
            g = g->next;
            } while ( e_start!=g) ; 
          }

          tf->eliminate();
          facemap.erase(tf);
        }
        if(f->exteriorEdge == te || f->exteriorEdge == e) f->exteriorEdge = tn;
        if(hole_c == 1 && is_hole(f->exteriorEdge)) f->exteriorEdge = n;
        if(hole_c == 2 && is_hole(f->exteriorEdge)) f->exteriorEdge = tn;

         // traverse all the deleted face and assigned the new face f
        
        tn->incidentFace = f;
        tp->incidentFace = f;
         HalfEdge* he = f->exteriorEdge ;
        const HalfEdge* e_start = he;
        do {
            he->incidentFace = f;
            he = he->next;
            } while ( e_start!=he) ; 

      }

      // we checked the hedge e, now check for other faces 
      //std::cout << "nin" <<  s.size()  << std::endl;  
      groupmap[e] =1; groupmap[te] =1;
      if(groupmap[n] == 0) s.push(n->twin);
      if(groupmap[p] == 0) s.push(p->twin);
      if(groupmap[tn] == 0) s.push(tn->twin);
      if(groupmap[tp] == 0) s.push(tp->twin);
      D.cleanup();
    }
    }
    D.cleanup();
  }
}
}
// 5.
int faces_in_mesh(std::unordered_map< Face*, int> &facemap, int i)
{
  int res = 0;
  for (auto & f : facemap )
  {
    if(f.second == i)
    {
      res ++;
    }
  }
  return res;
}


void exportCityJSON(DCEL & D, const char *file_out ,std::unordered_map< HalfEdge*, std::vector<int> > &hemap, std::unordered_map<int, Vertex*> &umap,
  std::unordered_map< Face*, int> &facemap) {

  std::ofstream myfile;
  std:: string  output = file_out;
  output = "../../" + output;
  myfile.open(output);
  std::cout << "Writing file:" << file_out << std::endl;
  myfile << "{";
  //init
  myfile <<"\"type\":\"CityJSON\",\"version\":\"1.0\",";
  //city objects
  myfile << "\"CityObjects\":{";
  int building = 1;
  for( const auto & e : D.infiniteFace()->holes)
  {
    myfile << "\"Building_" << building << "\":{";
    myfile << "\"geometry\":[{";
    myfile << "\"boundaries\":[";
    int last_f = 0 ;
    for( auto & f : facemap)
    {
      if(facemap[f.first] == building)
      {
        //boundaries start
        myfile << "[[";
        HalfEdge* e = f.first->exteriorEdge;
        const HalfEdge* e_start = e;
        int last = 0;
        do {
          int index;
          for( auto ver : umap)
          {
            if(ver.second == e->origin)
            {
              index = ver.first-1;
            }
          }
          if (e_start!=e->next) myfile << index << ",";
          else  myfile << index;
          last++;
          e = e->next;
        } while ( e_start!=e) ; 
        last_f ++;
        // print holes
        int i = 0;
        for (auto eh: f.first->holes)
        {
          std::cout << "Hole writing" << std::endl;
          //if the face has a hole, finish the exterior and start the interior
          myfile << "],[";
          const HalfEdge* e_starth = eh;
          int indexh = 0;
          do {
            for( auto ver : umap)
            {
              if(ver.second == eh->origin)
              {
                indexh = ver.first-1;
              }
            }
            if (e_starth!=eh->next) myfile << indexh << ",";
            else  myfile << indexh;
            eh = eh->next;
          } while ( e_starth!=eh) ; 
          // if the faces has 2 or more holes put a comma
          //if(i+1 <f.first->holes.size() ) myfile << ",";
          i++;
        }
        //finish the holes or the face
        myfile << "]";
      // finish the boundaries
      //std:: cout << last_f << "," << faces_in_mesh(facemap , building ) << std::endl;
      if(faces_in_mesh(facemap , building ) >  last_f) myfile << "],";
      else myfile << "]";
    }
    }
    //close boundaries
    myfile << "],";
    myfile << "\"lod\":2,\"type\":\"MultiSurface\"";
    //close geometry
    myfile   << "}],";
    myfile << "\"type\": \"Building\"";
    // close building
    myfile  << "}";
    if(building < D.infiniteFace()->holes.size()) myfile << "," ;
    building ++;
  }
  //close cityObjects
  myfile << "}, ";
  //vertices
  myfile << "\"vertices\":[";
  for (int i = 1; i < umap.size()+1; i++ ) {
  auto  v = umap[i];
  if(i == umap.size()){
    myfile << "["  << v->x << "," << v->y << "," <<v->z  << "]";
  }
  else myfile   <<  "[" <<v->x << "," <<v->y << "," <<v->z <<"]" <<",";
  }
  // close vertices
  myfile << "]";
  //close json
  myfile << "}";
  myfile.close();
  std::cout << "FINISHED Writing file: " << file_out << std::endl;
}


int main(int argc, const char * argv[]) {
  const char *file_in = "bk_soup.obj";
  const char *file_out = "bk.json";
  int count = 0;
  // Demonstrate how to use the DCEL to get you started (see function implementation below)
  // you can remove this from the final code
  //DemoDCEL();
 
  // create an empty DCEL
  DCEL D;
  //create an unordered map 
  std::unordered_map< HalfEdge*, std::vector<int> > hemap;
  std::unordered_map<int, Vertex*> umap;

  //  create unrdered map

  // 1. read the triangle soup from the OBJ input file and convert it to the DCEL,
  importOBJ(D, file_in,  hemap, umap);
  // 2. group the triangles into meshes,
  std::unordered_map< Face*, int> facemap;
  for( auto & f : facemap)
  {
    facemap.insert({f.first,0});
  }
   groupTriangles(D, hemap, facemap);
  // 3. determine the correct orientation for each mesh and ensure all its triangles 
  //    are consistent with this correct orientation (ie. all the triangle normals 
  //    are pointing outwards).
  std::unordered_map< Face*, int> fmap;
  for( auto & f : facemap)
  {
    fmap.insert({f.first,0});
  }
   OrientationOfMap(D, fmap, facemap);
    std::unordered_map< HalfEdge*, int> groupmap;
 
  for( auto & f : groupmap)
  {
    groupmap.insert({f.first,0});
  }
  // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
  mergeCoPlanarFaces( D,facemap, groupmap);
  // 5. write the meshes with their faces to a valid CityJSON output file.
  exportCityJSON(D, file_out, hemap, umap, facemap);
  return 0;
}


