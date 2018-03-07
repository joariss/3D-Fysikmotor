/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pkg3dstuff;

/**
 *
 * @author joari
 */
public class TestStuff {
    static double[][] sphereNodes = {{0,0,0}};
    public static void main(String[] args) {
        int radius = 10;
        
        for (int i = 0; i < radius*2; i++) {
            
            double r2 = Math.sqrt(8*i*radius) / 2 - Math.sqrt(4*Math.pow(i,2)) / 2;
            for (double v = 0; v < r2*Math.PI/(Math.PI / 180); v+= 1) {
                AllClasses.addDoubleArray(sphereNodes, new double[] {Math.cos(v/(Math.PI / 180)) * r2, i, Math.sin(v/(Math.PI / 180)) * r2});
                
            }
        }
        System.out.println(sphereNodes);
    }
}
