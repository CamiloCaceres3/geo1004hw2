#include <iostream>
#include <list>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stack>

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
  std::unordered_map< HalfEdge*, int> meshmap;
  for( auto & e : hemap)
  {
    meshmap.insert({e.first,0});
  }
  int build = 1;
  bool list_em = false;
  while(list_em == false)
  {
    HalfEdge* c;
    std::unordered_map<HalfEdge*, int>::iterator it = meshmap.begin();
    while (it != meshmap.end()) {
        if (it->second == 0) {
            c = it->first;
            it == meshmap.end();
        } 
        it++;
    }

    std::stack<HalfEdge*> s;
    s.push(c);
    D.infiniteFace()->holes.push_back(c);
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
    for(auto c: meshmap)
    {
      if(c.second ==0)
      {
        ir ++;
      }
    }
    if(ir == 0)
    {
      list_em = true;
    }
  }

std:: cout << D.infiniteFace()->holes.size() << std::endl;


}
// 3.
void orientMeshes(DCEL & D) {
  // to do
}
// 4.
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

bool coplanar(Vertex* &v1, Vertex* &v2, Vertex* &v3, Vertex* &v4)
{
  std::vector<double> V1 {v2->x - v1->x, v2->y - v1->y , v2->z -v1->z };
  std::vector<double> V2 {v4->x - v1->x, v4->y - v1->y , v4->z -v1->z };
  std::vector<double> V3 {v3->x - v1->x, v3->y - v1->y , v3->z -v1->z };

  std::vector<double> v  = cross(V1, V2);
  double d =  dot(v, V3);

  if(abs(d) < 0.01)
  {
    return true;
  }
  else return false;
}


void mergeCoPlanarFaces(DCEL & D,   std::unordered_map< Face*, int> &facemap) {
  // to do
  for(  auto hedge: D.infiniteFace()->holes)
  {
    HalfEdge* e = hedge;
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

    //check if the are coplanar
    bool cop = coplanar(v1,v2,v3,v4);
    std::cout << cop << *v1 << *v2 << *v3 << *v4 <<std::endl;

    if(cop == true )
    {
      //if(e == hedge ) D.infiniteFace()->holes.push_back(n);
      te->incidentFace->eliminate();
      e -> eliminate();
      te -> eliminate ();
      facemap.erase(te->incidentFace);

      n->prev = tp;
      p->next = tn;
      tn->prev = p;
      tp->next = n;
      tn->incidentFace = f;
      tp->incidentFace = f;
      f->exteriorEdge = n;
    }

    DCELElement* I  = D.findInValid();
    //Vertex* vv = HalfEdge I->destination;
    std::cout <<  *e << *te << *n << *p << *tn << *tp <<std::endl;
    printDCEL(D);
    D.cleanup();
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
        myfile << "]]";
        if(faces_in_mesh(facemap , building ) >  last_f) myfile << ",";
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
    myfile   << *v ;
  }
  else myfile   << *v <<",";
  }
  // close vertices
  myfile << "]";
  //close json
  myfile << "}";
  myfile.close();
  std::cout << "FINISHED Writing file: " << file_out << std::endl;
}


int main(int argc, const char * argv[]) {
  const char *file_in = "polygonal_hole.obj";
  const char *file_out = "bk.json";

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
  
  // 4. merge adjacent triangles that are co-planar into larger polygonal faces.
  mergeCoPlanarFaces( D,facemap);
  // 5. write the meshes with their faces to a valid CityJSON output file.
  exportCityJSON(D, file_out, hemap, umap, facemap);
  return 0;
}


void printDCEL(DCEL & D) {

  // Quick check if there is an invalid element
  auto element = D.findInValid();
  if ( element == nullptr ) {
    // Beware that a 'valid' DCEL here only means there are no dangling links and no elimated elements.
    // There could still be problems like links that point to the wrong element.
    std::cout << "DCEL is valid\n";
  } else {
    std::cout << "DCEL is NOT valid ---> ";
    std::cout << *element << "\n";
  }

  // iterate all elements of the DCEL and print the info for each element
  const auto & vertices = D.vertices();
  const auto & halfEdges = D.halfEdges();
  const auto & faces = D.faces();
  std::cout << "DCEL has:\n";
  std::cout << " " << vertices.size() << " vertices:\n";
  for ( const auto & v : vertices ) {
    std::cout << "  * " << *v << "\n";
  }
  std::cout << " " << halfEdges.size() << " half-edges:\n";
  for ( const auto & e : halfEdges ) {
    std::cout << "  * " << *e << "\n";
  }
  std::cout << " " << faces.size() << " faces:\n";
  for ( const auto & f : faces ) {
    std::cout << "  * " << *f << "\n";
  }

}


void DemoDCEL() {

  std::cout << "/// STEP 1 Creating empty DCEL...\n";
  DCEL D;
  printDCEL(D);

  /*

  v2 (0,1,0)
   o
   |\
   | \
   |  \
   o---o v1 (1,0,0)
  v0
  (0,0,0)

  We will construct the DCEL of a single triangle 
  in the plane z=0 (as shown above).

  This will require:
    3 vertices
    6 halfedges (2 for each edge)
    1 face

  */
  std::cout << "\n/// STEP 2 Adding triangle vertices...\n";
  Vertex* v0 = D.createVertex(0,0,0);
  Vertex* v1 = D.createVertex(1,0,0);
  Vertex* v2 = D.createVertex(0,1,0);
  printDCEL(D);

  std::cout << "\n/// STEP 3 Adding triangle half-edges...\n";
  HalfEdge* e0 = D.createHalfEdge();
  HalfEdge* e1 = D.createHalfEdge();
  HalfEdge* e2 = D.createHalfEdge();
  HalfEdge* e3 = D.createHalfEdge();
  HalfEdge* e4 = D.createHalfEdge();
  HalfEdge* e5 = D.createHalfEdge();
  printDCEL(D);

  std::cout << "\n/// STEP 4 Adding triangle face...\n";
  Face* f0 = D.createFace();
  printDCEL(D);

  std::cout << "\n/// STEP 5 Setting links...\n";
  e0->origin = v0;
  e0->destination = v1;
  e0->twin = e3;
  e0->next = e1;
  e0->prev = e2;
  e0->incidentFace = f0;

  e3->origin = v1;
  e3->destination = v0;
  e3->twin = e0;
  e3->next = e5;
  e3->prev = e4;

  /* 
  If a half-edge is incident to 'open space' (ie not an actual face with an exterior boundary), 
  we use the infiniteFace which is predifined in the DCEL class
  */
  e3->incidentFace = D.infiniteFace();

  e1->origin = v1;
  e1->destination = v2;
  e1->twin = e4;
  e1->next = e2;
  e1->prev = e0;
  e1->incidentFace = f0;

  e4->origin = v2;
  e4->destination = v1;
  e4->twin = e1;
  e4->next = e3;
  e4->prev = e5;
  e4->incidentFace = D.infiniteFace();

  e2->origin = v2;
  e2->destination = v0;
  e2->twin = e5;
  e2->next = e0;
  e2->prev = e1;
  e2->incidentFace = f0;

  e5->origin = v0;
  e5->destination = v2;
  e5->twin = e2;
  e5->next = e4;
  e5->prev = e3;
  e5->incidentFace = D.infiniteFace();

  f0->exteriorEdge = e0;
  printDCEL(D);


  std::cout << "\n/// STEP 6 Traversing exterior vertices of f0...\n";
  /* 
  if all is well in the DCEL, following a chain of half-edges (ie keep going to e.next)
  should lead us back the the half-edge where we started.
  */
  HalfEdge* e = f0->exteriorEdge;
  const HalfEdge* e_start = e;
  do {
    std::cout << " -> " << *e->origin << "\n";
    e = e->next;
  } while ( e_start!=e) ; // we stop the loop when e_start==e (ie. we are back where we started)
  
  
  std::cout << "\n/// STEP 7 eliminating v0...\n";
  v0->eliminate();
  printDCEL(D);
  
  /* 
  We just eliminated v0. At the same time we know there are elements that still 
  pointers to v0 (ie the edges e0, e2, e3, e5). This means we can NOT call D.cleanup()!
  If you do this anyways, the program may crash. 
  
  Eg. if you uncomment the following there could be a crash/stall of the program.
  */
  // D.cleanup(); // this will remove v0 from memory (because we just eliminated v0 and the cleanup() function simply removes all the eliminated elements)
  // std::cout << *v0 << "\n"; // we try to access that memory, but v0 is gone -> undefined behaviour 
  // std::cout << *e0->origin << "\n"; // this equivalent to the previous line (both point to the same memory address)


  std::cout << "\n/// STEP 8 eliminating all the remaining DCEL elements\n";
  for ( const auto & v : D.vertices() ) {
    v->eliminate();
  }
  for ( const auto & e : D.halfEdges() ) {
    e->eliminate();
  }
  for ( const auto & f : D.faces() ) {
    f->eliminate();
  }
  printDCEL(D);

  std::cout << "\n/// STEP 9 cleaning up the DCEL\n";
  D.cleanup();
  printDCEL(D);

}
