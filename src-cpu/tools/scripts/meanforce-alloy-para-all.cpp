#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <chrono>


int user_ntype;
double user_rcut;
int user_nbin;

// mean force for alloy and only for dump style file


using namespace std;
class atom{
    public:
        double x,y,z;
        double fx,fy,fz;
        double xs,ys,zs;
        int type;
};         //定义类

class meanf{
        
   public:
   void val(int m_atomi, int m_atomj, int user_nbin);
   int mfCalu(int user_nbin);
   int atomi;
   int atomj;
   double* mf;
   double* count;
   double* mf_ave;
   double* gr_ave;
   double* mfInt;
};         //定义类

void read_para(char* argv[]);






int main(int argc, char* argv[]) {               //主函数请求两个值，一个值给argc，另一个存储所有的输入参数，argv[0]存程序的名称
   auto start = chrono::high_resolution_clock::now();
    if( argc != 3)
    {
       cout << "exe para atomfile" << endl;
       return 0;
    }
    cout << "Reading para.in" << endl;
    read_para(argv);
    int nstep=0;  
    extern double user_rcut;                  //cut-off
    extern int user_ntype;                    //初始化原子种类数目
    extern int user_nbin;
    int npfcf=user_ntype*(user_ntype+1)/2;
    double box[3];                              //box是一个数组
    double bl, bh;                              //box size
    int natom;                                  //总原子数                            //nbin数目
    double dr=user_rcut/user_nbin;                        //dr
    atom* atoms;                                //创建atom类的一个指针
    meanf* mfs;
    atom* atoms2;
    ifstream fin(argv[2]);                      //从文件中定义实例:fin
    string ss;                                  //定义字符串变量
    double ii;                                  //
    int timestep;
    cout << "******" << endl;                   //
    cout << "Computing mean force for "   << argv[2] << endl; //
    cout << "rcut: " << user_rcut << "; dr: "<< dr << "; user_nbin: " <<  user_nbin   <<  endl;      
    
   mfs=new meanf [npfcf];
   for (int i = 0; i < user_ntype; i++)
   {
      for (int j = i; j < user_ntype; j++)
      {
         int index=(2*i*user_ntype-i*(i-1))/2+j-i;
         cout << "ppcf for pair "   <<j << i   << " " << index <<" is made "  << endl; 
         mfs[index].val(i+1,j+1,user_nbin);
      }
   }
   cout << "making ppcfs "  << endl; 

    while (1) 
    {               //每一个snaps都进行计算

       getline(fin,ss);
       fin >> timestep; getline(fin,ss);getline(fin,ss);         //读取timestep
       fin >> natom ; getline(fin,ss);          //读取atom但是没跳行 
      //  cout << " MD step "<< timestep <<endl;         
       getline(fin,ss);                         //跳过一行    
       for (int k=0; k<3; k++) {                //读取bound信息
          fin >> bl >> bh ; getline(fin,ss);    //   
          box[k]=bh-bl;                         //
       }
       getline(fin,ss);                         //dump style
       
       // get atomic position and force, compute MF
       atoms = new atom [natom];                //new返回第一个元素的地址   

       
       for (int k=0; k<natom; k++) {            //创建完了所有的atoms实例
          fin >> ii >> atoms[k].type 
              >> atoms[k].xs >> atoms[k].ys >> atoms[k].zs
              >> atoms[k].fx >> atoms[k].fy >> atoms[k].fz ;
          getline(fin,ss);
       }
      //  cout << " MD step "<< atoms[-1].xs <<endl;
       for (int i=0; i<natom-1; i++)
          for (int j=i+1; j<natom; j++) {
             double dx=atoms[i].xs-atoms[j].xs;
             double dy=atoms[i].ys-atoms[j].ys;
             double dz=atoms[i].zs-atoms[j].zs;
             if(dx>0.5) dx=1-dx;
             if(dy>0.5) dy=1-dy;
             if(dz>0.5) dz=1-dz;
             if(dx<-0.5) dx=1+dx;
             if(dy<-0.5) dy=1+dy;
             if(dz<-0.5) dz=1+dz;
             dx=dx*box[0]; dy=dy*box[1]; dz=dz*box[2];
             double dis=sqrt( dx*dx + dy*dy +dz*dz ); 
             if(dis>user_rcut) continue;
             //if(dis< 1.0) {cout << "Warning: distance less than 1 A! " <<  endl;}
             double dfx=atoms[i].fx-atoms[j].fx;   //i粒子合力-j粒子合力 所以当i粒子为a,j粒子为b时，对平均力的贡献是phi_ab，如果粒子i为b，粒子j为a时，对平均力的贡献是phi_ba，按照结果，这两者不会有区别？
             double dfy=atoms[i].fy-atoms[j].fy;
             double dfz=atoms[i].fz-atoms[j].fz;
             
             // Eqn. 8
             double fr=0.5*(dfx*dx+dfy*dy+dfz*dz)/dis ;

             int m=int(dis/dr);         //感觉也可以把所有前面的count加一个值
             int index_j = max(atoms[i].type,atoms[j].type)-1;   int index_i = min(atoms[i].type,atoms[j].type) -1 ;
             int typeCount=(2*index_i*user_ntype-index_i*(index_i-1))/2+index_j-index_i;
             //cout  << "atom pair "  << index_j  << "  " << index_i <<  "belong to " <<typeCount  << endl;
            //  if (typeCount>=npfcf)
            //  { 
            //    cout<<"------"<<endl;
            //    cout  << "Number of the types is not matches."   << endl;
            //    cout  << "******"   << endl;
            //    return 0;
            //  }
             mfs[typeCount].mf[m] +=fr ;
             mfs[typeCount].count[m] +=1;
         
          }
      //  cout  << "nstep: " << nstep <<endl;
       delete [] atoms;
       nstep+=1;         //完全存疑
       if (fin.peek() == EOF ) 
         {
            cout  << "file is all readed, exit " <<endl;
               // 如果文件已经被读取完毕，退出循环
               break;
         }
         
         // cout   <<" step " << nstep << " atoms " << natom <<" MD step "<< timestep <<endl;
    }
    fin.close();
    cout<<"------"<<endl;

    int nmin=99999, rmin=0;


   cout<<"------"<<endl;
   ofstream fout1("mean-force-alloy.dat");
   fout1 <<"# r ";
   for (int i = 0; i <npfcf ; i++)
      {
         rmin = mfs[i].mfCalu(user_nbin);
         cout  << "meanforce type: "<< mfs[i].atomj  << mfs[i].atomi 
         << " bin with least pop. at r=" << dr*(rmin+0.5) 
         << ", Npop=" << mfs[i].count[rmin] << endl;
         fout1 << "mean-force-"  <<  mfs[i].atomj   << mfs[i].atomi   <<" " ;
      }
   fout1  << endl;


   for (int k = 0; k < user_nbin; k++)
   {
      fout1 << dr*(k+0.5);
      for (int i = 0; i <npfcf ; i++)
      {
         fout1 << " " << mfs[i].mf_ave[k];
      }
      fout1 << endl;
   }

   fout1.close();
   delete [] mfs;
   cout << "******" << endl;
   cout  << "nstep of all steps " << nstep <<endl;
   auto end = std::chrono::high_resolution_clock::now();

    // 计算程序运行时间
   std::chrono::duration<double> duration = end - start;

   std::cout << "Program execution time: " << duration.count() << " seconds\n";

   return 0;
   
}
void read_para(char* argv[]){
   ifstream fp(argv[1]);
   string ss,filename;
   string job;
   while (fp>>ss)
   {
      if (ss=="NTYPE")
         fp >> user_ntype;
      if (ss=="NBIN")
         fp >> user_nbin;
      if (ss=="RCUT")
         fp >> user_rcut;
   }
   fp.close();
}
void meanf::val(int m_atomi,int m_atomj,int user_nbin){
   atomi=m_atomi;
   atomj=m_atomj;
   mf=new double[user_nbin];
   count=new double[user_nbin];
   mf_ave=new double[user_nbin];
   gr_ave=new double[user_nbin];
   for (int k=0; k<user_nbin; k++) {
         mf[k]=0;
         mf_ave[k]=0;
         count[k]=0;
         gr_ave[k]=0;
   } 
}
int meanf::mfCalu(int user_nbin){

   int nmin=99999, rmin=0;
      double pi=3.14159265359;
      for (int k=0; k<user_nbin; k++) {
         mf_ave[k]=mf[k]/count[k];
         if(count[k]!=0 && count[k]<nmin) {
            nmin=count[k];
            rmin=k;
         }
      }        
   return rmin;

}