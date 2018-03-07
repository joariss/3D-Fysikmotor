package pkg3dstuff;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import javax.swing.JPanel;

public class MainPanel extends JPanel implements Runnable, MouseListener, MouseMotionListener, KeyListener{

BufferedImage canvas = new BufferedImage(1920, 1080, BufferedImage.TYPE_INT_ARGB);




public MainPanel() {
    
    new Thread(this).start();
    addMouseListener(this);
    addMouseMotionListener(this);
    addKeyListener(this);
    
}


public void run(){
    while (true) {
        requestFocus();
        repaint();
        try {
            Thread.sleep(1);
        } catch (Exception e) {
        }
    }
}

Color[] originalColors = {new Color(100,100,100)};
static double[][] originalNodes = {{0,0,0}};
static int[][] originalEdges = {};
static double[][] originalSphereNodes = {{0,0,0}};
Color[] colors = originalColors;
double[][] SphereNodes = originalSphereNodes;
double[][] nodes = originalNodes;

int[][] edges = originalEdges;

double mouseX, mouseY, pMouseX, pMouseY;
double mouseXMovedDist;
double mouseYMovedDist;

double mouseMovedXTotal;
double mouseMovedYTotal;

double[][] sphereDistanceMoved = {
    {0,0,0},
    {0,0,0},
    {0,0,0}, 
    {0,0,0},
};

double[][] screenDistanceMoved = {{0,0,0}};

double[][] sphereSpeed = {
    {0,0,0},
    {0,0,0},
    {0,0,0}, 
    {0,0,0}, 
};

int a = 0;
static boolean forwardPressed = false;
static boolean rightPressed = false;
static boolean backwardPressed = false;
static boolean leftPressed = false;
static boolean upPressed = false;
static boolean downPressed = false;


void draw() {
    Graphics g = canvas.getGraphics();
    g.setColor(Color.WHITE);
    g.fillRect(0,0,1920,1080);
    
    //sphereInfo[i][0] = position       (double[]);
    //sphereinfo[i][1] = radius         (double);
    //sphereInfo[i][2] = acceleration       (double[]);
    //sphereInfo[i][3] = speed          (double[]);
    //sphereInfo[i][4] = distanceMoved; (double[]);
    //sphereInfo[i][5] = mass           (double);
    //////sphereInfo[i][5] = rotation           (double);
    
    double[][][] sphereInfo = {
    //    0                 1      2          3                  4                5
        {{ 0, -250, 0}, {50}, {0,0,0}, sphereSpeed[0], sphereDistanceMoved[0], {1}},  
        {{ 0, -100, 0}, {50}, {0,0,0}, sphereSpeed[1], sphereDistanceMoved[1], {1}}, 
        {{ 0, 0, 0}, {50}, {0,0,0}, sphereSpeed[2], sphereDistanceMoved[2], {4}}, 
        //{{ 100, -200, 500}, {50}, {0,0,0}, sphereSpeed[3], sphereDistanceMoved[3], {1}}, 
        
    };
        Object.StuckSpring(sphereInfo[0], new double[]{-100,-400,0}, 1, 1, 100, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[1], sphereInfo[2], 0.1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[0], sphereInfo[1], 10, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.Spring(sphereInfo[1], sphereInfo[0], 1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{100,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{500,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); AllClasses.doASDFSPACESHIFTPressed(sphereDistanceMoved, screenDistanceMoved);

        Object.StuckSpring(sphereInfo[0], new double[]{200,-400,0}, 1, 1, 100, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[1], sphereInfo[2], 0.1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[0], sphereInfo[1], 10, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.Spring(sphereInfo[1], sphereInfo[0], 1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{100,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{500,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); AllClasses.doASDFSPACESHIFTPressed(sphereDistanceMoved, screenDistanceMoved);

        Object.Spring(sphereInfo[0], sphereInfo[1], 1, 1, 100, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[1], sphereInfo[2], 0.1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[0], sphereInfo[1], 10, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.Spring(sphereInfo[1], sphereInfo[0], 1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{100,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{500,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); AllClasses.doASDFSPACESHIFTPressed(sphereDistanceMoved, screenDistanceMoved);
        Object.Spring(sphereInfo[2], sphereInfo[1], 1, 1, 100, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[1], sphereInfo[2], 0.1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); // Object.Spring(sphereInfo[0], sphereInfo[1], 10, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.Spring(sphereInfo[1], sphereInfo[0], 1, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{100,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); Object.StuckSpring(sphereInfo[0], new double[]{500,-300,100}, 0.3, 1, 1, mouseMovedXTotal, mouseMovedYTotal , g); AllClasses.doASDFSPACESHIFTPressed(sphereDistanceMoved, screenDistanceMoved);

 for (int i = 0; i < sphereInfo.length; i++) { 
        for (int j = 0; j < i; j++) { 
        if(AllClasses.dist3D(new double[] {sphereInfo[i][0][0] + sphereInfo[i][4][0], sphereInfo[i][0][1] + sphereInfo[i][4][1], sphereInfo[i][0][2] + sphereInfo[i][4][2]}, new double[] {sphereInfo[j][0][0] + sphereInfo[j][4][0], sphereInfo[j][0][1] + sphereInfo[j][4][1], sphereInfo[j][0][2] + sphereInfo[j][4][2]}) < (sphereInfo[i][1][0] + sphereInfo[j][1][0])){
            
            sphereDistanceMoved[j] = AllClasses.checkSphereColision12(sphereInfo[j], sphereInfo[i]);
            sphereDistanceMoved[i] = AllClasses.checkSphereColision22(sphereInfo[j], sphereInfo[i]);
            sphereSpeed[j] = AllClasses.checkSphereColision11(sphereInfo[j], sphereInfo[i]);
            sphereSpeed[i] = AllClasses.checkSphereColision21(sphereInfo[j], sphereInfo[i]);
            
        } 
    } 
} 

    //System.out.println(sphereInfo[0][3][0]);
    
    AllClasses.groundColission(sphereInfo, screenDistanceMoved);
    
    
    
    
    for (int i = 0; i < sphereInfo.length; i++) {
        AllClasses.addGravity(sphereInfo[i]);
        Object.Sphere(new double[][] {sphereInfo[i][0]}, sphereInfo[i][4] ,Object.CreateSphere((int)sphereInfo[i][1][0], SphereNodes), sphereInfo[i][2], mouseMovedXTotal, mouseMovedYTotal , colors, g); 

        AllClasses.accelerateSphere(sphereInfo[i]);
        
        AllClasses.moveSphere(sphereInfo[i]);
    }

    
}


public void paint(Graphics g) {
    draw();
    g.drawImage(canvas, 0, 0, this);
}

@Override
public void mouseClicked(MouseEvent e) {
    
}

@Override
public void mousePressed(MouseEvent e) {
    mouseX = e.getX();
    mouseY = e.getY();
    pMouseX = mouseX;
    pMouseY = mouseX;
}

@Override
public void mouseReleased(MouseEvent e) {
    
}

@Override
public void mouseEntered(MouseEvent e) {
}

@Override
public void mouseExited(MouseEvent e) {
}

    @Override
    public void mouseDragged(MouseEvent e) {

        pMouseX = mouseX;
        pMouseY = mouseY;
        mouseX = e.getX();
        mouseY = e.getY();




        mouseXMovedDist = AllClasses.mouseYMovedDist(e.getY(), pMouseY);
        mouseYMovedDist = AllClasses.mouseXMovedDist(e.getX(), pMouseX);

        mouseMovedXTotal += mouseXMovedDist;
        mouseMovedYTotal += mouseYMovedDist;
    // AllClasses.rotateX(mouseXMovedDist / 360 * 3 , nodes); // AllClasses.rotateY(mouseYMovedDist / 360 * 3 , nodes); }
    }
    @Override
    public void mouseMoved(MouseEvent e) {
    }

    @Override
    public void keyTyped(KeyEvent e) {
        if(e.getKeyCode() == KeyEvent.VK_UP){

        }
    }

    @Override
    public void keyPressed(KeyEvent e) {
       AllClasses.checkKeyPressedForMovement(e);
    }

    @Override
    public void keyReleased(KeyEvent e) {
        AllClasses.checkKeyReleasedForMovement(e);
    }
}