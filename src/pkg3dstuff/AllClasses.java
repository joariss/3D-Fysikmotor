/*

To change this license header, choose License Headers in Project Properties.
To change this template file, choose Tools | Templates
and open the template in the editor. */ package pkg3dstuff;
import java.awt.Color; import java.awt.Graphics; import java.awt.List; import java.awt.event.KeyEvent; import java.util.ArrayList; import java.util.Arrays; import static pkg3dstuff.MainPanel.forwardPressed;
import java.lang.Math.*;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static java.lang.Math.tan;
/** *

@author joagra458 */ public class AllClasses { 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    static Color[] addColors(Color[] colors, Color addedColor){ Color[] newColorArr= Arrays.copyOf(colors, colors.length+1); newColorArr[colors.length]= addedColor; colors= newColorArr; return colors; }


















static int[][] addEdges(int[][] edges, int[] addedEdge){ if(addedEdge[0] > 1){ int[][] newEdgeArr= Arrays.copyOf(edges, edges.length+1); newEdgeArr[edges.length]= addedEdge; edges= newEdgeArr; } return edges; }

static double[][] addDoubleArray(double[][] nodes, double[] addedNode){
    double[][] newNodeArr= Arrays.copyOf(nodes, nodes.length+1);
    newNodeArr[nodes.length]= addedNode; nodes= newNodeArr;
    return nodes;
};

static double[][][] addSecondDoubleArray(double[][][] nodes, double[][] addedNode){ double[][][] newNodeArr= Arrays.copyOf(nodes, nodes.length+1); newNodeArr[nodes.length]= addedNode; nodes= newNodeArr; return nodes; };

static void drawEdges(double x, double y, double z, int[][] edges, double[][] nodes, Graphics g, Color[] color){ int fov = 1920; for (int i=0; i<edges.length; i++) { int n0 = edges[i][0]; int n1 = edges[i][1]; double[] node0 = nodes[n0]; double[] node1 = nodes[n1]; g.setColor(color[i]);

     int p00 = (int)(fov * (node0[0] + x) / (fov + node0[2] + z));
     int p01 = (int)(fov * (node0[1] + y) / (fov + node0[2] + z));
     
     int p10 = (int)(fov * (node1[0] + x) / (fov + node1[2] + z));
     int p11 = (int)(fov * (node1[1] + y) / (fov + node1[2] + z));
     
     //if(node0[2] + z > 0 && node1[2] + z > 0){
         g.drawLine(p00 + Main.width / 2, p01 + Main.height / 2, p10 + Main.width / 2, p11 + Main.height / 2);
     //}
     
 }
}

static void drawSprings(double[] position1, double[] position2, Graphics g){ int fov = 1920; int p1x = (int)(fov * (position1[0]) / (fov + position1[2])); int p1y = (int)(fov * (position1[1]) / (fov + position1[2]));

 int p2x = (int)(fov * (position2[0]) / (fov + position2[2]));
 int p2y = (int)(fov * (position2[1]) / (fov + position2[2]));
 
 g.setColor(Color.BLACK);
 
 g.drawLine(p1x + Main.width / 2, p1y + Main.height / 2, p2x + Main.width / 2, p2y + Main.height / 2);
}

static void drawFaces(int[][] faces, double[][] nodes, Graphics g, Color[] color){ int fov = 1920; for(int i = 0; i < faces.length; i ++){ if(faces[i].length == 3){ int n0 = faces[i][0]; int n1 = faces[i][1]; int n2 = faces[i][2];

         double[] node0 = nodes[n0];
         double[] node1 = nodes[n1];
         double[] node2 = nodes[n2];
         
         int p00 = (int)(fov * node0[0] / (fov + node0[2]));
         int p01 = (int)(fov * node0[1] / (fov + node0[2]));
         
         int p10 = (int)(fov * node1[0] / (fov + node1[2]));
         int p11 = (int)(fov * node1[1] / (fov + node1[2]));
         
         int p20 = (int)(fov * node2[0] / (fov + node2[2]));
         int p21 = (int)(fov * node2[1] / (fov + node2[2]));
         
         g.drawLine(p00 + 700, p01 + 700, p10 + 700, p11 + 700);
         g.drawLine(p10 + 700, p11 + 700, p20 + 700, p21 + 700);
         g.drawLine(p20 + 700, p21 + 700, p00 + 700, p01 + 700);
         
         int a = (int)(Math.sqrt(Math.pow(p00 - p10, 2) + Math.pow(p01 + p11, 2))+ 0.5) ;
         int b = (int)(Math.sqrt(Math.pow(p10 - p20, 2) + Math.pow(p11 + p21, 2))+ 0.5) ;
         int c = (int)(Math.sqrt(Math.pow(p20 - p00, 2) + Math.pow(p21 + p01, 2))+ 0.5) ;
         
         
         
         
     }
     
     if(faces[i].length == 4){
         int n0 = faces[i][0];
         int n1 = faces[i][1];
         int n2 = faces[i][2];
         int n3 = faces[i][3];
         
         double[] node0 = nodes[n0];
         double[] node1 = nodes[n1];
         double[] node2 = nodes[n2];
         double[] node3 = nodes[n3];
         
         int p00 = (int)(fov * node0[0] / (fov + node0[2]));
         int p01 = (int)(fov * node0[1] / (fov + node0[2]));
         
         int p10 = (int)(fov * node1[0] / (fov + node1[2]));
         int p11 = (int)(fov * node1[1] / (fov + node1[2]));
         
         int p20 = (int)(fov * node2[0] / (fov + node2[2]));
         int p21 = (int)(fov * node2[1] / (fov + node2[2]));
         
         int p30 = (int)(fov * node3[0] / (fov + node3[2]));
         int p31 = (int)(fov * node3[1] / (fov + node3[2]));
         
         
         for (int j = 0; j < (int)p11 - p31; j++) {
             g.drawLine((int)(p00 + 0.5) + Main.width / 2,     (int)(p01 + 0.5) - j + Main.height / 2 , (int)(p10 + 0.5) + Main.width / 2,    (int)(p11 + 0.5) - j + Main.height / 2);
         }
         
     }
 }
}

static double mouseXMovedDist(double mouseX, double pMouseX){ double mouseXDistanceMoved = pMouseX - mouseX; return mouseXDistanceMoved; }

static double mouseYMovedDist(double mouseY, double pMouseY){ double mouseYDistanceMoved = pMouseY - mouseY; return mouseYDistanceMoved; }

public static double random(double min, double max) { return min + Math.random() * (max - min); }

static void rotateX(double theta, double[][] nodes){ double sinTheta = Math.sin(theta); double cosTheta = Math.cos(theta); for (int i=0; i<nodes.length; i++) { double [] node = nodes[i]; double y = node[1]; double z = node[2]; node[1] = y * cosTheta - z * sinTheta; node[2] = z * cosTheta + y * sinTheta; } }

static void rotateY(double theta, double[][] nodes){ double sinTheta = Math.sin(theta); double cosTheta = Math.cos(theta);

 for (int i=0; i<nodes.length; i++) {
     double [] node = nodes[i];
     double x = node[0];
     double z = node[2];
     node[0] = x * cosTheta - z * sinTheta;
     node[2] = z * cosTheta + x * sinTheta;
 }
}

static void rotateZ(double theta, double[][] nodes){ double sinTheta = Math.sin(theta); double cosTheta = Math.cos(theta);

 for (int n=0; n<nodes.length; n++) {
     double[] node = nodes[n];
     double x = node[0];
     double y = node[1];
     node[0] = x * cosTheta - y * sinTheta;
     node[1] = y * cosTheta + x * sinTheta;
 }
}

static void drawSphere(double[] spherePosition, double[][] sphereNodes, Graphics g){ int fov = 1920; g.setColor(Color.black); for (int i = 1; i < sphereNodes.length; i++) { double[] node0 = sphereNodes[i];

         int x = (int)(fov * (node0[0] + spherePosition[0]) / (fov + node0[2] + spherePosition[2]));
         int y = (int)(fov * (node0[1] + spherePosition[1]) / (fov + node0[2] + spherePosition[2]));
     
     if(i < sphereNodes.length - 20){
         double[] node1 = sphereNodes[i + 20];
         double[] node2 = sphereNodes[i + 1];
         if((i % 20) + 1== 1 && i != 1){
            node2 = sphereNodes[i - 19];
         }
         int x2 = (int)(fov * (node1[0] + spherePosition[0]) / (fov + node1[2] + spherePosition[2]));
         int y2 = (int)(fov * (node1[1] + spherePosition[1]) / (fov + node1[2] + spherePosition[2]));
         
         int x3 = (int)(fov * (node2[0] + spherePosition[0]) / (fov + node2[2] + spherePosition[2]));
         int y3 = (int)(fov * (node2[1] + spherePosition[1]) / (fov + node2[2] + spherePosition[2]));

         g.drawLine(x + Main.width / 2, y + Main.height / 2, x2 + Main.width / 2, y2 + Main.height / 2 );
         g.drawLine(x + Main.width / 2, y + Main.height / 2, x3 + Main.width / 2, y3 + Main.height / 2 );

     
     }

 }
}

static double dist1D(double a, double b){ return Math.abs(a - b); }

static double dist3D(double[] sphere1Position, double[] sphere2Position){
    return Math.sqrt(Math.pow((sphere1Position[0] - sphere2Position[0]), 2) + Math.pow((sphere1Position[1] - sphere2Position[1]), 2) + Math.pow((sphere1Position[2] - sphere2Position[2]), 2)); }





static void moveCamera(int moveX, int moveY, int moveZ, double[][]sphereDistanceMoved){ for (int i = 0; i < sphereDistanceMoved.length; i++) { sphereDistanceMoved[i][0] += moveX; sphereDistanceMoved[i][1] += moveY; sphereDistanceMoved[i][2] += moveZ; } }

//static double[] sphereCollision(double[][] sphere1, double[][]sphere2){
//    double r12 = sphere1[1][0] + sphere2[1][0];
//    
//    double m21 = sphere2[5][0] / sphere1[5][0];
//    
//    double x1 = sphere1[0][0] + sphere1[4][0];
//    double y1 = sphere1[0][0] + sphere1[4][0];
//    double z1 = sphere1[0][0] + sphere1[4][0];
//    
//    double x21 = x1 - (sphere1[0][0] + sphere1[4][0]);
//    double y21 = y1 - (sphere1[0][1] + sphere1[4][1]);
//    double z21 = z1 - (sphere1[0][2] + sphere1[4][2]);
//    
//    double vx21 = sphere2[3][0] - sphere1[3][0];
//    double vy21 = sphere2[3][1] - sphere1[3][1];
//    double vz21 = sphere2[3][2] - sphere1[3][2];
//    
//    double d = Math.sqrt(x21*x21 + y21*y21 + z21*z21);
//    double v = Math.sqrt(vx21*vx21 + vy21*vy21 + vz21*vz21);
//    
//    double error = 0;
//    
//    
//    double x2 = x21;
//    double y2 = y21;
//    double z2 = z21;
//    
//    double vx1 = -vx21;
//    double vy1 = -vy21;
//    double vz1 = -vz21;
//    double phi2 = 0;
//    
//    double vx2 = sphere2[3][0];
//    double vy2 = sphere2[3][1];
//    double vz2 = sphere2[3][2];
//    
//    double theta2 = Math.acos(z2 / d);
//    if(x2 == 0 && y2 == 0){
//        phi2 = 0;
//    }else{
//        phi2 = Math.atan2(y2, x2);
//    }
//    
//    double st = Math.sin(theta2);
//    double ct = Math.cos(theta2);
//    double sp = Math.sin(phi2);
//    double cp = Math.cos(phi2);
//    
//    double vx1r = ct*cp*vx1+ct*sp*vy1-st*vz1;
//    double vy1r = cp*vy1-sp*vx1;
//    double vz1r = st*cp*vx1+st*sp*vy1+ct*vz1;
//    
//    double fvz1r = vz1r/v;
//    
//    if(fvz1r > 1){
//        fvz1r = 1;
//    }else if(fvz1r < -1){
//        fvz1r = -1;
//    }
//    double phiv = 0;
//    double thetav = Math.acos(fvz1r);
//    if(vx1r == 0 && vy1r == 0){
//        phiv = 0;
//    }else{
//        phiv = Math.atan2(vy1r, vx1r);
//    }
//    
//    double dr = d*Math.sin(thetav)/r12;
//    
//    double alpha = Math.asin(-dr);
//    double beta = phiv;
//    double sbeta=Math.sin(beta);
//    double cbeta=Math.cos(beta);
//    
//    double t = (d*Math.cos(thetav) - r12*Math.sqrt(1-dr*dr))/v;
//    
//    x2 = x2+vx2*t+x1;
//    y2 = y2+vy2*t+y1;
//    z2 = z2+vz2*t+z1;
//    
//    x1 =(vx1+vx2)*t +x1;
//    y1 =(vy1+vy2)*t +y1;
//    z1 =(vz1+vz2)*t +z1;
//    
//    double a = Math.tan(thetav + alpha);
//    
//    double dvz2 = 2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*1+m21);
//    
//    double vz2r = dvz2;
//    double vx2r = a*cbeta*dvz2;
//    double vy2r = a*sbeta*dvz2;
//    
//    vz1r = vz1r-m21*vz2r;
//    vx1r = vx1r-m21*vx2r;
//    vy1r = vy1r-m21*vy2r;
//    
//   vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
//   vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
//   vz1=ct*vz1r-st*vx1r               +vz2;
//   vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
//   vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
//   vz2=ct*vz2r-st*vx2r               +vz2;
//   
//   double[] b = new double[]{vx1, vy1, vz1, vx2, vy2, vz2};
//   if(d<r12){
//        b = new double[]{sphere1[3][0], sphere1[3][1], sphere1[3][2], sphere2[3][0], sphere2[3][1], sphere2[3][2]};
//    }
//   if(v == 0){
//        b = new double[]{sphere1[3][0], sphere1[3][1], sphere1[3][2], sphere2[3][0], sphere2[3][1], sphere2[3][2]};
//    }
//   return b;
//    
//}


static double[] sphereCollision(double[][] sphere1, double[][]sphere2){
    double x1 = sphere1[0][0] + sphere1[4][0];
    double y1 = sphere1[0][1] + sphere1[4][1];
    double z1 = sphere1[0][2] + sphere1[4][2];
    double x2= sphere2[0][0] + sphere2[4][0];
    double y2= sphere2[0][1] + sphere2[4][1];
    double z2= sphere2[0][2] + sphere2[4][2];
                                     
    double vx1 = sphere1[3][0];
    double vy1 = sphere1[3][1];
    double vz1 = sphere1[3][2];
                                             
    double vx2 = sphere2[3][0];
    double vy2 = sphere2[3][1];
    double vz2 = sphere2[3][2];
    int error;


       double  pi,r12,m21,d,v,theta2,phi2,st,ct,sp,cp,vx1r,vy1r,vz1r,fvz1r,
	           thetav,phiv,dr,alpha,beta,sbeta,cbeta,dc,sqs,t,a,dvz2,
			   vx2r,vy2r,vz2r,x21,y21,z21,vx21,vy21,vz21,vx_cm,vy_cm,vz_cm;

//     **** initialize some variables ****
       double r1 = sphere1[1][0];
       double r2 = sphere2[1][0];
       double m1 = sphere1[5][0];
       double m2 = sphere2[5][0];
       pi=acos(-1.0E0);
       error=0;
       r12=r1+r2;
       m21=m2/m1;
       x21=x2-x1;
       y21=y2-y1;
       z21=z2-z1;
       vx21=vx2-vx1;
       vy21=vy2-vy1;
       vz21=vz2-vz1;
       
       vx_cm = (m1*vx1+m2*vx2)/(m1+m2) ;
       vy_cm = (m1*vy1+m2*vy2)/(m1+m2) ;
       vz_cm = (m1*vz1+m2*vz2)/(m1+m2) ;  

	   
//     **** calculate relative distance and relative speed ***
       d=sqrt(x21*x21 +y21*y21 +z21*z21);
       v=sqrt(vx21*vx21 +vy21*vy21 +vz21*vz21);
       
//     **** return if distance between balls smaller than sum of radii ****
       
//     **** return if relative speed = 0 ****
       
       

//     **** shift coordinate system so that ball 1 is at the origin ***
       x2=x21;
       y2=y21;
       z2=z21;
       
//     **** boost coordinate system so that ball 2 is resting ***
       vx1=-vx21;
       vy1=-vy21;
       vz1=-vz21;

//     **** find the polar coordinates of the location of ball 2 ***
       theta2=acos(z2/d);
       if (x2==0 && y2==0) {phi2=0;} else {phi2=atan2(y2,x2);}
       st=sin(theta2);
       ct=cos(theta2);
       sp=sin(phi2);
       cp=cos(phi2);


//     **** express the velocity vector of ball 1 in a rotated coordinate
//          system where ball 2 lies on the z-axis ******
       vx1r=ct*cp*vx1+ct*sp*vy1-st*vz1;
       vy1r=cp*vy1-sp*vx1;
       vz1r=st*cp*vx1+st*sp*vy1+ct*vz1;
       fvz1r = vz1r/v ;
       if (fvz1r>1) {fvz1r=1;}   // fix for possible rounding errors
          else if (fvz1r<-1) {fvz1r=-1;} 
       thetav=acos(fvz1r);
       if (vx1r==0 && vy1r==0) {phiv=0;} else {phiv=atan2(vy1r,vx1r);}

        						
//     **** calculate the normalized impact parameter ***
       dr=d*sin(thetav)/r12;


//     **** return old positions and velocities if balls do not collide ***
//       if (thetav>pi/2 || fabs(dr)>1) {
//           x2=x2+x1;
//           y2=y2+y1;
//           z2=z2+z1;
//           vx1=vx1+vx2;
//           vy1=vy1+vy2;
//           vz1=vz1+vz2;
//           error=1;
//           return;
//        }
       
//     **** calculate impact angles if balls do collide ***
       alpha=asin(-dr);
       beta=phiv;
       sbeta=sin(beta);
       cbeta=cos(beta);
        
       
//     **** calculate time to collision ***
       t=(d*cos(thetav) -r12*sqrt(1-dr*dr))/v - 1;

     
//     **** update positions and reverse the coordinate shift ***
//       x2=x2+vx2*t +x1 - sphere2[0][0];
//       y2=y2+vy2*t +y1 - sphere2[0][1];
//       z2=z2+vz2*t +z1 - sphere2[0][2];
//       x1=(vx1+vx2)*t +x1 - sphere1[0][0];
//       y1=(vy1+vy2)*t +y1 - sphere1[0][1];
//       z1=(vz1+vz2)*t +z1 - sphere1[0][2];
       
       x2=x2+vx2*t +x1 - sphere2[0][0];
       y2=y2+vy2*t +y1 - sphere2[0][1];
       z2=z2+vz2*t +z1 - sphere2[0][2];
       x1=(vx1+vx2)*t +x1 - sphere1[0][0];
       y1=(vy1+vy2)*t +y1 - sphere1[0][1];
       z1=(vz1+vz2)*t +z1 - sphere1[0][2];
       
//       x2= MainPanel.oldPositions[1][0];
//       y2= MainPanel.oldPositions[1][1];
//       z2= MainPanel.oldPositions[1][2];
//       x1= MainPanel.oldPositions[0][0];
//       y1= MainPanel.oldPositions[0][1];
//       z1= MainPanel.oldPositions[0][2];
        
 
       
//  ***  update velocities ***

       a=tan(thetav+alpha);

       dvz2=2*(vz1r+a*(cbeta*vx1r+sbeta*vy1r))/((1+a*a)*(1+m21));
       
       vz2r=dvz2;
       vx2r=a*cbeta*dvz2;
       vy2r=a*sbeta*dvz2;
       vz1r=vz1r-m21*vz2r;
       vx1r=vx1r-m21*vx2r;
       vy1r=vy1r-m21*vy2r;

       
//     **** rotate the velocity vectors back and add the initial velocity
//           vector of ball 2 to retrieve the original coordinate system ****
                     
       vx1=ct*cp*vx1r-sp*vy1r+st*cp*vz1r +vx2;
       vy1=ct*sp*vx1r+cp*vy1r+st*sp*vz1r +vy2;
       vz1=ct*vz1r-st*vx1r               +vz2;
       vx2=ct*cp*vx2r-sp*vy2r+st*cp*vz2r +vx2;
       vy2=ct*sp*vx2r+cp*vy2r+st*sp*vz2r +vy2;
       vz2=ct*vz2r-st*vx2r               +vz2;
        

//     ***  velocity correction for inelastic collisions ***

//       vx1=(vx1-vx_cm)*R + vx_cm;
//       vy1=(vy1-vy_cm)*R + vy_cm;
//       vz1=(vz1-vz_cm)*R + vz_cm;
//       vx2=(vx2-vx_cm)*R + vx_cm;
//       vy2=(vy2-vy_cm)*R + vy_cm;
//       vz2=(vz2-vz_cm)*R + vz_cm;  
        double[] b = new double[]{vx1, vy1, vz1, vx2, vy2, vz2, x1, y1, z1, x2, y2, z2};
        x21 = x2-x1;
        y21 = y2-x1;
        z21 = z2-x1;
        d = sqrt(x21*x21+y21*y21+z21*z21);
       if (d<r12) {
        error=2;
            //System.out.println(error);
            //b = new double[]{sphere1[3][0], sphere1[3][1], sphere1[3][2], sphere2[3][0], sphere2[3][1], sphere2[3][2], sphere1[4][0], sphere1[4][1], sphere1[4][2], sphere2[4][0], sphere2[4][1], sphere2[4][2]};
       }
       if (v==0) {
        error=1;
       //b = new double[]{sphere1[3][0], sphere1[3][1], sphere1[3][2], sphere2[3][0], sphere2[3][1], sphere2[3][2], sphere1[4][0], sphere1[4][1], sphere1[4][2], sphere2[4][0], sphere2[4][1], sphere2[4][2]};
       }
       return b;
}

static double[] checkSphereColision11(double[][] sphere1, double[][] sphere2){ 
    double[] a = {};
        a = sphereCollision(sphere1, sphere2); 
     
    double[] b = new double[]{a[0], a[1], a[2]};
    
    return b;
}
static double[] checkSphereColision12(double[][] sphere1, double[][] sphere2){ 
    double[] a = {};
        a = sphereCollision(sphere1, sphere2); 
     
    double[] b = new double[]{a[6], a[7], a[8]};
    
    return b;
}
static double[] checkSphereColision21(double[][] sphere1, double[][] sphere2){ 
    double[] a = {};
        a = sphereCollision(sphere1, sphere2); 
    
    double[] b = new double[]{a[3], a[4], a[5]};
    return b;
}
static double[] checkSphereColision22(double[][] sphere1, double[][] sphere2){ 
    double[] a = {};
        a = sphereCollision(sphere1, sphere2); 
    
    double[] b = new double[]{a[9], a[10], a[11]};
    return b;
}

static void addSpeedToSphere(double[][][] sphereInfo){

}

static void doASDFSPACESHIFTPressed(double[][] sphereDistanceMoved, double[][] screenDistanceMoved){ if(MainPanel.forwardPressed){ AllClasses.moveCamera(0,0,-1,sphereDistanceMoved); AllClasses.moveCamera(0,0,-1,screenDistanceMoved); } if(MainPanel.rightPressed){ AllClasses.moveCamera(-1,0,0,sphereDistanceMoved); AllClasses.moveCamera(-1,0,0,screenDistanceMoved); } if(MainPanel.backwardPressed){ AllClasses.moveCamera(0,0,1,sphereDistanceMoved); AllClasses.moveCamera(0,0,1,screenDistanceMoved); } if(MainPanel.leftPressed){ AllClasses.moveCamera(1,0,0,sphereDistanceMoved); AllClasses.moveCamera(1,0,0,screenDistanceMoved); } if(MainPanel.upPressed){ AllClasses.moveCamera(0,1,0,sphereDistanceMoved); AllClasses.moveCamera(0,1,0,screenDistanceMoved); } if(MainPanel.downPressed){ AllClasses.moveCamera(0,-1,0,sphereDistanceMoved); AllClasses.moveCamera(0,-1,0,screenDistanceMoved); } }

static void groundColission(double[][][] sphereInfo, double[][] screenDistanceMoved){ for (int i = 0; i < sphereInfo.length; i++) { if(dist1D(sphereInfo[i][0][1] + sphereInfo[i][4][1], screenDistanceMoved[0][1] + 300) < sphereInfo[i][1][0]){ sphereInfo[i][3][1] = -sphereInfo[i][3][1]; } } }

static void addGravity(double[][] sphereInfo){ addForceToSPhere(sphereInfo, new double[]{0,9.81/1000 * sphereInfo[5][0],0}); }

static void checkKeyPressedForMovement(KeyEvent e){ if(e.getKeyCode() == KeyEvent.VK_W){ MainPanel.forwardPressed = true; } if(e.getKeyCode() == KeyEvent.VK_D){ MainPanel.rightPressed = true; } if(e.getKeyCode() == KeyEvent.VK_S){ MainPanel.backwardPressed = true; } if(e.getKeyCode() == KeyEvent.VK_A){ MainPanel.leftPressed = true; } if(e.getKeyCode() == KeyEvent.VK_SPACE){ MainPanel.upPressed = true; } if(e.getKeyCode() == KeyEvent.VK_SHIFT){ MainPanel.downPressed = true; } }

static void checkKeyReleasedForMovement(KeyEvent e){ if(e.getKeyCode() == KeyEvent.VK_W){ MainPanel.forwardPressed = false; } if(e.getKeyCode() == KeyEvent.VK_D){ MainPanel.rightPressed = false; } if(e.getKeyCode() == KeyEvent.VK_S){ MainPanel.backwardPressed = false; } if(e.getKeyCode() == KeyEvent.VK_A){ MainPanel.leftPressed = false; } if(e.getKeyCode() == KeyEvent.VK_SPACE){ MainPanel.upPressed = false; } if(e.getKeyCode() == KeyEvent.VK_SHIFT){ MainPanel.downPressed = false; } }

static void addForceToSPhere(double[][] sphereInfo, double[] addedForce){ sphereInfo[2][0] += addedForce[0] / sphereInfo[5][0]; sphereInfo[2][1] += addedForce[1] / sphereInfo[5][0]; sphereInfo[2][2] += addedForce[2] / sphereInfo[5][0]; }

static void accelerateSphere(double[][] sphereInfo){ sphereInfo[3][0] += sphereInfo[2][0]; sphereInfo[3][1] += sphereInfo[2][1]; sphereInfo[3][2] += sphereInfo[2][2]; }

static void moveSphere(double[][] sphereInfo){ sphereInfo[4][0] += sphereInfo[3][0]; sphereInfo[4][1] += sphereInfo[3][1]; sphereInfo[4][2] += sphereInfo[3][2]; }

static void addSpringForce(double[][] sphere1, double[][] sphere2, double springConstant, double friction, double springLength){ double yxv = Math.atan2(sphere2[0][1] - sphere1[0][1], sphere2[0][0]-sphere1[0][0]); double zyv = Math.atan2(sphere2[0][2] - sphere1[0][2], sphere2[0][1]-sphere1[0][1]); double distanceBetweenSpheresX = (sphere2[0][0]) - (sphere1[0][0]); double distanceBetweenSpheresY = (sphere2[0][1]) - (sphere1[0][1]); double distanceBetweenSpheresZ = (sphere2[0][2]) - (sphere1[0][2]);

 addForceToSPhere(sphere1, new double[]{springConstant *  distanceBetweenSpheresX / 1000, springConstant *  distanceBetweenSpheresY / 1000, springConstant *  distanceBetweenSpheresZ / 1000});
 addForceToSPhere(sphere2, new double[]{springConstant * -distanceBetweenSpheresX / 1000, springConstant * -distanceBetweenSpheresY / 1000, springConstant * -distanceBetweenSpheresZ / 1000});
}

static void addStuckSpringForce(double[][] sphere1, double[] StuckPosition, double springConstant, double friction, double springLength){

 double yxv = Math.atan2(sphere1[0][1] - StuckPosition[1], StuckPosition[0]-sphere1[0][0]);
 double zyv = Math.atan2(sphere1[0][2] - StuckPosition[2], StuckPosition[1]-sphere1[0][1]);
 double distanceBetweenSpheresX = (StuckPosition[0]) - (sphere1[0][0]);
 double distanceBetweenSpheresY = (StuckPosition[1]) - (sphere1[0][1]);
 double distanceBetweenSpheresZ = (StuckPosition[2]) - (sphere1[0][2]);
// double distanceBetweenSpheresX = (StuckPosition[0]) - (sphere1[0][0]); // double distanceBetweenSpheresY = (StuckPosition[1]) - (sphere1[0][1]); // double distanceBetweenSpheresZ = (StuckPosition[2]) - (sphere1[0][2]);

    addForceToSPhere(sphere1, new double[]{springConstant *  distanceBetweenSpheresX / 1000, springConstant *  distanceBetweenSpheresY / 1000, springConstant *  distanceBetweenSpheresZ / 1000});
}
}