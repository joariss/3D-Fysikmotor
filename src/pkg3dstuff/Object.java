/*

To change this license header, choose License Headers in Project Properties.
To change this template file, choose Tools | Templates
and open the template in the editor. */ package pkg3dstuff;
import java.awt.Color; import java.awt.Graphics; import java.awt.List; import java.util.ArrayList;

/** *

@author joagra458 */ public class Object { static void Object(double[][] objectPosition, double[][] nodes, int[][] edges, int[][] faces, double mouseXMovedDist,double mouseYMovedDist, Color[] colors, Graphics g){
// AllClasses.rotateX(mouseXMovedDist / 360 * 3 , objectPosition); // AllClasses.rotateY(mouseYMovedDist / 360 * 3 , objectPosition); AllClasses.rotateX(mouseXMovedDist / 360 * 3 , nodes); AllClasses.rotateY(mouseYMovedDist / 360 * 3 , nodes); for(int i = 0; i < edges.length; i ++){ colors = AllClasses.addColors(colors, new Color(0,0,0)); }

    AllClasses.drawEdges(objectPosition[0][0], objectPosition[0][1], objectPosition[0][2], edges, nodes, g, colors);
}

static double[][] CreateSphere(int radius, double[][] sphereNodes){
    int asd = 10;
    for (double i = 0; i <= radius*2; i+= radius/asd) {
        double r2 = Math.sqrt(8*radius*i - 4*Math.pow(i, 2)) / 2;
            for (double v = 0; v < 2*Math.PI; v+= Math.PI / asd) {
                sphereNodes = AllClasses.addDoubleArray(sphereNodes, new double[] {(Math.cos(v))  * r2 , i - radius, (Math.sin(v)) * r2});
            }
    }
    return sphereNodes;
}


static void Sphere(double[][] spherePosition, double[] distanceMoved, double[][] sphereNodes, double[] objectRotation, double mouseXMovedDist,double mouseYMovedDist, Color[] colors, Graphics g){

    spherePosition[0][0] += distanceMoved[0];
    spherePosition[0][1] += distanceMoved[1];
    spherePosition[0][2] += distanceMoved[2];

    AllClasses.rotateX(mouseXMovedDist / 360 * 3, spherePosition);
    AllClasses.rotateY(mouseYMovedDist / 360 * 3, spherePosition);

    AllClasses.rotateX((objectRotation[0] + mouseXMovedDist) / 360 * 3, sphereNodes);
    AllClasses.rotateY((objectRotation[1] + mouseYMovedDist) / 360 * 3, sphereNodes);
    AllClasses.rotateZ(objectRotation[2], sphereNodes);

    AllClasses.drawSphere(spherePosition[0], sphereNodes, g);
    
    spherePosition[0][0] -= distanceMoved[0];
    spherePosition[0][1] -= distanceMoved[1];
    spherePosition[0][2] -= distanceMoved[2];
}

static void Spring(double[][] sphere1, double[][] sphere2, double springConstant, double friction, double springLength, double mouseXMovedDist,double mouseYMovedDist, Graphics g){
    
    sphere1[0][0] += sphere1[4][0];
    sphere1[0][1] += sphere1[4][1];
    sphere1[0][2] += sphere1[4][2];
    
    sphere2[0][0] += sphere2[4][0] ;
    sphere2[0][1] += sphere2[4][1];
    sphere2[0][2] += sphere2[4][2];
    
    double[][] Position = {new double[]{sphere1[0][0], sphere1[0][1], sphere1[0][2]}, new double[]{sphere2[0][0], sphere2[0][1], sphere2[0][2]}};
    
    AllClasses.addSpringForce(sphere1, sphere2, springConstant, friction, springLength);
    
    AllClasses.rotateX(mouseXMovedDist / 360 * 3, Position);
    AllClasses.rotateY(mouseYMovedDist / 360 * 3, Position);
    
    
    AllClasses.drawSprings(Position[0], Position[1], g);
    
    
    sphere2[0][0] -= sphere2[4][0];
    sphere2[0][1] -= sphere2[4][1];
    sphere2[0][2] -= sphere2[4][2];
    
    sphere1[0][0] -= sphere1[4][0];
    sphere1[0][1] -= sphere1[4][1];
    sphere1[0][2] -= sphere1[4][2];
}

static void StuckSpring(double[][] sphere, double[] springPosition, double springConstant, double friction, double springLength, double mouseXMovedDist,double mouseYMovedDist, Graphics g){
    sphere[0][0] += sphere[4][0] ;
    sphere[0][1] += sphere[4][1];
    sphere[0][2] += sphere[4][2];
    
    double[][] Position = {new double[]{sphere[0][0], sphere[0][1], sphere[0][2]}, springPosition};
    
    AllClasses.addStuckSpringForce(sphere, springPosition, springConstant, friction, springLength);
    
    AllClasses.rotateX(mouseXMovedDist / 360 * 3, Position);
    AllClasses.rotateY(mouseYMovedDist / 360 * 3, Position);
    
    
    AllClasses.drawSprings(Position[0], springPosition, g);
    
    
    sphere[0][0] -= sphere[4][0] ;
    sphere[0][1] -= sphere[4][1];
    sphere[0][2] -= sphere[4][2];
}
}

