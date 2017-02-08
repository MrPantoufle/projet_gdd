Edge S; // the surface, aka initial handle on an edge
float distance ;  // for camera
int[] vColor = {100, 0, 0};
int[] nColor = {0, 100, 0};
int vertexSize = 10;
int maxFaceDegree = 8; // maximum number of edges in a face (suboptimal, used in ImportPLY)
float scale;

float f = 100;

void setup() {
  size(500, 500, P3D);
  distance = width;
  S = importPLY("dodecahedron.ply");
  //S = pyramid();
  scale = width/10;
}


void draw() {
  background(50);
  // horizontal angle phi = PI*(2*mouseX/width -1)
  // vertical angle theta = (PI/2)*(2*mouseY/height/2 - 1)
  camera(width/2-distance*sin(TWO_PI*mouseX/width)*sin(PI*mouseY/height), 
    height/2+distance*cos(PI*mouseY/height), 
    -distance*cos(TWO_PI*mouseX/width)*sin(PI*mouseY/height), 
    width/2, height/2, 0, 
    0, 1, 0);



  ArrayList<Face> list = faceTraversal(S.face);
  for (Face currentFace : list) {
    currentFace.draw(250, 137);
  }

  ArrayList<Vertex> vlist = vertexTraversal(S.first);
  for (Vertex currentVertex : vlist) {
    ArrayList<Vertex> neigh = currentVertex.neighbors();
    int red = 0;
    Vector h = new Vector();
    for(Vertex v : neigh) {
       Vector tmp = new Vector(v.xpos, v.ypos, v.zpos);
       h.add(tmp);
    }
    red = (int)(h.norm()/(scale*10)*255);
    int[] colors = {red, 128, 0};
    currentVertex.drawVertex(colors);
  }
  
  ArrayList<Vector> flow = badCurvatureFlow(vlist);
  S.shift(flow);
}

void mouseWheel(MouseEvent event) {  // for zooming in and out
  float e = event.getCount();
  distance += e;
}

ArrayList<Vector> badCurvatureFlow(ArrayList<Vertex> vertices) {
  Vector centre = gravity(vertices);
  
  ArrayList<Vector> shifts = new ArrayList<Vector>();
  for(Vertex currentVer : vertices) {
    Vector currentShift = new Vector();
    ArrayList<Vertex> neigh = currentVer.neighbors();
    for(Vertex v : neigh) {
       Vector tmp = new Vector(v.xpos, v.ypos, v.zpos);
       //tmp.normalize();
       currentShift.add(tmp);
    }
    //currentShift.normalize();
    //currentShift.mult(0.1);
    shifts.add(currentShift);
  }
  
  for(Vertex currentVer : vertices) {
      currentVer.xpos -= centre.x;
      currentVer.ypos -= centre.y;
      currentVer.zpos -= centre.z;
  }

  
  return shifts;
}

Vector gravity(ArrayList<Vertex> vertices) {
   Vector centre = new Vector();
   for(Vertex v : vertices) {
      centre.x += v.xpos;
      centre.y += v.ypos;
      centre.z += v.zpos;
   }
   centre.mult(1.0/(float)vertices.size());
  return centre;
}


///////////////// CLASSES /////////////////

class Vertex {
  float xpos;
  float ypos;
  float zpos;
  Edge edgefrom;   // one of the oriented edges emanating from the vertex
  int status;    // for listing purposes

  Vertex(float x, float y, float z) {  
    xpos = x + width/2;
    ypos = y + height/2;
    zpos = z;
    status = 0;
  }

  void drawVertex(int[] fillColor) {
    //noStroke();
    stroke(fillColor[0], fillColor[1], fillColor[2]);
    lights();
    pushMatrix();
    translate(xpos, ypos, zpos);
    sphere(vertexSize);
    popMatrix();
  }

  ArrayList<Vertex> neighbors() {  // ordered list of adjacent vertices
    ArrayList<Vertex> nlist = new ArrayList<Vertex>();
    Edge current = this.edgefrom;
    Vertex start = current.last;   // adjacent vertex at the end of edgefrom 
    do {
      nlist.add(current.last);
      current = current.previous.opposite;
    } while (current.last != start);
    return nlist;
  }
}


class Edge {  // oriented
  Vertex first;
  Vertex last;
  Edge opposite;
  Edge next;
  Edge previous;
  Face face;
  int status;    // for listing purposes

  Edge(Vertex theFirst, Vertex theLast) {
    first = theFirst;
    last = theLast;
    status = 0;
  }

  void draw() {
    stroke(255);
    line(first.xpos, first.ypos, first.zpos, last.xpos, last.ypos, last.zpos);
  }
  
  void shift(ArrayList<Vector> shifts) {
    ArrayList<Vertex> list = vertexTraversal(this.first);
    // one should check that the size of both list is the same ...
    //Vector v,s;
    float biggestNorm = 0;
    // We return the biggest shift vector, to
    // have an adaptative flow.
    for (int i = 0; i < list.size(); i++) {
      Vector s = shifts.get(i);
      biggestNorm = max(biggestNorm, s.norm());
    }
    biggestNorm = width/(f*biggestNorm);
    for(int i = 0; i < list.size(); i++) {
      Vector s = shifts.get(i);
      s.mult(biggestNorm);
      list.get(i).xpos += s.x;
      list.get(i).ypos += s.y;
      list.get(i).zpos += s.z;
    }
    //return biggestNorm;
  }
}


class Face { //oriented
  Edge inedge;  // one (oriented) edge
  int degree;  // face degree, i.e. number of edges
  int status;    // for listing purposes
  String name;

  Face(Edge e, int d) {
    inedge = e;
    degree = d;
    status = 0;
  }

  void draw(int theStroke, int theFill) {
    if (theStroke == -1) {
      noStroke();
    } else {
      stroke(theStroke);
    }
    if (theFill == -1) {
      noFill();
    } else {
      fill(theFill);
    }
    Edge e = inedge;
    beginShape();
    for (int i=0; i <= degree; i++) {
      vertex(e.last.xpos, e.last.ypos, e.last.zpos);
      e = e.next;
    }
    endShape();
  }

  void drawEdges(int strokeColor) {
    draw(strokeColor, -1);
  }

  void drawFace(int fillColor) {
    draw(-1, fillColor);
  }

  ArrayList<Edge> boundingEdges() {  // returns the list of boundary edges
    Edge e = inedge;
    ArrayList<Edge> list = new ArrayList<Edge>();
    for (int i=0; i < degree; i++) {
      list.add(e);
      e = e.next;
    }
    return(list);
  }

  ArrayList<Face> incidentFaces() {  // returns the list of incident faces
    ArrayList<Face> faces = new ArrayList<Face>();
    Edge e = inedge;
    for (int i=0; i < degree; i++) {
      faces.add(e.opposite.face);
      e = e.next;
    }
    return(faces);
  }
}

class Vector {
  float x = 0.0;
  float y = 0.0;
  float z = 0.0;

  Vector() {
  }

  Vector(float xi, float yi, float zi) {
    x = xi;
    y = yi;
    z = zi;
  }

  void add(Vector v) {
    x += v.x;
    y += v.y;
    z += v.z;
  }
  void sub(Vector v) {
    x-= v.x;
    y-= v.y;
    z-= v.z;
  }

  void mult(float c) {
    x = x*c;
    y = y*c;
    z = z*c;
  }

  float det(Vector v1, Vector v2) {
     return x*(v1.y*v2.z - v1.z*v2.y) - y*(v1.x*v2.z - v1.z*v2.x) + z*(v1.x*v2.y - v1.y*v2.x);
  }

  float scal(Vector v) {
    return (x * v.x + y *v.y + z*v.z);
  }

  float norm() {
    return sqrt(normsq());
  }
  float normsq() {
    return x*x+y*y+z*z;
  }
  void normalize() {
    float mag = this.normsq();
    x = x/mag;
    y = y/mag;
    z = z/mag;
  }
}


class Surface {
  // not as half-edge convention
  Edge initialEdge;

  Surface(Edge e) {
    initialEdge = e;
  }

  void draw() {
    // draws all faces by enumeration
  }
  
  void shift(ArrayList<Vector> shifts) {
    ArrayList<Vertex> list = vertexTraversal(this.initialEdge.first);
    // one should check that the size of both list is the same ...
    //Vector v,s;
    for (int i = 0; i < list.size(); i++) {
      Vector s = shifts.get(i);
      list.get(i).xpos += s.x;
      list.get(i).ypos += s.y;
      list.get(i).zpos += s.z;
    }
  }
}


///////////  ROUTINES   //////////////////

Edge importPLY(String filename) {
  String[] lines = loadStrings(filename);
  int numVertices = 0;
  int numFaces = 0;
  boolean end_header = false;
  String plyTest = "ply";
  if (!lines[0].equals(plyTest)) exit();
  float scalingFactor = width/5;

  int i = 0;
  while (!end_header) {
    String[] keywords = split(lines[i], ' ');
    if (keywords[0].equals("element")) {
      if (keywords[1].equals("vertex")) {
        numVertices = int(keywords[2]);
      } else if (keywords[1].equals("face")) {
        numFaces = int(keywords[2]);
      }
    } else if (keywords[0].equals("end_header")) {
      end_header = true;
    }
    i++;
  }
  println("v=", numVertices, " f=", numFaces);

  Vertex[] v = new Vertex[numVertices];
  Edge[] e = new Edge[2*maxFaceDegree*numFaces];  // two directions!
  Face[] f = new Face[numFaces];

  // Vertex' 3D coordinates
  for (int j = 0; j < numVertices; j++) {
    String[] keywords = split(lines[i], ' ');
    v[j] = new Vertex(scalingFactor*float(keywords[0]), scalingFactor*float(keywords[1]), scalingFactor*float(keywords[2]));
    i++;    // increase line number
  }

  // list faces and obtain combinatorial structure
  // store the edges in a temporary Edge-valued incidence matrix
  // so that edges are indexed by two vertex indices 

  // start by creating a complete graph, to garantee valid pointers
  Edge[][] edgeIncidence = new Edge[numVertices][numVertices];
  for (int j = 0; j < numVertices; j++) {
    for (int k = 0; k < numVertices; k++) {
      edgeIncidence[j][k] = new Edge(v[j], v[k]);
    }
  }
  for (int j = 0; j < numFaces; j++) {  // list all faces 
    String[] keywords = split(lines[i], ' ');
    int degree = int(keywords[0]);
    // set the degree and status of the face with first edge as inedge
    f[j] = new Face(null, degree);  
    f[j].name = "f"+char(j);

    for (int k = 0; k < degree; k++) { // list through the edges
      int startV = int(keywords[1+k]);
      int endV = int(keywords[1+(k + 1)%degree]);

      edgeIncidence[startV][endV].next = edgeIncidence[endV][int(keywords[1+(k+2)%degree])];
      edgeIncidence[startV][endV].previous = edgeIncidence[int(keywords[1+(k+degree-1)%degree])][startV];
      edgeIncidence[startV][endV].face = f[j];
      // add vertex edgefrom (vertices may and will be overwritten)
      v[startV].edgefrom = edgeIncidence[startV][endV];
    }
    f[j].inedge = edgeIncidence[int(keywords[1])][int(keywords[2])];
    i++; // increase line number
  }
  // fill opposite edge, if it exists (i.e. has valid face info)
  // this is necessary for surfaces with boundary
  for (int j = 0; j < numVertices; j++) {
    for (int k = j+1; k < numVertices; k++) {
      if (edgeIncidence[k][j].face != null) {  // opposite edge really exists
        edgeIncidence[j][k].opposite = edgeIncidence[k][j];
        edgeIncidence[k][j].opposite = edgeIncidence[j][k];
      }
    }
  }

  // transform t edge 2D array into an edge 1D array
  // only edges with "face" field are valid
  int currentVIndex = 0;  // index in array e
  for (i = 0; i < numVertices; i++) {
    for (int j = 0; j < numVertices; j++) {
      if (edgeIncidence[i][j].face != null) {
        e[currentVIndex] = edgeIncidence[i][j];
        currentVIndex++;
      }
    }
  }

  return(f[0].inedge);
}



ArrayList<Face> faceTraversal(Face initialFace) {  // lists all faces starting from initialFace
  // assumes all face statuses are = 0
  // uses Dijkstra-like algorithm
  ArrayList<Face> finalList = new ArrayList<Face>();
  ArrayList<Face> discoveredList = new ArrayList<Face>();  // list of discovered faces
  ArrayList<Face> neighborsList = new ArrayList<Face>();
  Face currentFace;
  //Edge currentEdge = initialFace.inedge;
  discoveredList.add(initialFace);
  while (discoveredList.size() > 0) {
    currentFace = discoveredList.get(0);
    discoveredList.remove(0);
    neighborsList = currentFace.incidentFaces();
    for (Face f : neighborsList) {
      if (f.status == 0) {   // new face
        f.status = 1;   // switch to discovered
        discoveredList.add(f);
      }
    }
    finalList.add(currentFace);
    currentFace.status = 2;    // switch to marked
  }
  // clean all faces once done
  for (Face f : finalList) {
    f.status = 0;
  }
  return(finalList);
}


ArrayList<Vertex> vertexTraversal(Vertex initialVertex) {  // lists all vertices starting from initialVertex
  // assumes all vertex statuses are = 0
  // uses Dijkstra-like algorithm
  ArrayList<Vertex> finalList = new ArrayList<Vertex>();
  ArrayList<Vertex> discoveredList = new ArrayList<Vertex>();  // list of discovered vertices
  ArrayList<Vertex> neighborsList = new ArrayList<Vertex>();
  Vertex currentVertex;
  discoveredList.add(initialVertex);
  while (discoveredList.size() > 0) {
    currentVertex = discoveredList.get(0);
    discoveredList.remove(0);
    neighborsList = currentVertex.neighbors();
    for (Vertex f : neighborsList) {
      if (f.status == 0) {   // new Vertex
        f.status = 1;   // switch to discovered
        discoveredList.add(f);
      }
    }
    finalList.add(currentVertex);
    currentVertex.status = 2;    // switch to marked
  }
  // clean all faces once done
  for (Vertex f : finalList) {
    f.status = 0;
  }
  return(finalList);
}