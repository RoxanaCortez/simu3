#include <fstream>

void obtenerDatos(istream &file,int nlines,int n,int mode,item* item_list){
    string line;
    file >> line;

    if(nlines==DOUBLELINE)
      file >> line;

   for(int i=0;i<n;i++){
        switch(mode){
        case INT_FLOAT:
            int e; float r;
            file >> e >> r;
            item_list[i].setIntFloat(e,r);
            break;
        case INT_INT_INT:
            int e1,e2,e3;
            file >> e1 >> e2 >> e3;
            item_list[i].setIntIntInt(e1,e2,e3);
            break;
        }
    }
}

void leerMallayCondiciones(mesh &m){
    char filename[10];
    ifstream file;
    float e,a,f;
    int nnodes,neltos,ndirich,nneu;

  do{
     cout << "Ingrese el nombre del archivo que contiene los datos de la malla: ";
     cin >> filename;
     file.open(filename);
    }while(!file);

    file >> e >> a >> f;
    file >> nnodes >> neltos >> ndirich >> nneu;

    m.setParameters(e,a,f);
    m.setSizes(nnodes,neltos,ndirich,nneu);
    m.createData();

    obtenerDatos(file,SINGLELINE,nnodes,INT_FLOAT,m.getNodes());
    obtenerDatos(file,DOUBLELINE,neltos,INT_INT_INT,m.getElements());
    obtenerDatos(file,DOUBLELINE,ndirich,INT_FLOAT,m.getDirichlet());
    obtenerDatos(file,DOUBLELINE,nneu,INT_FLOAT,m.getNeumann());

    file.close();
}


