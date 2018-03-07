package pkg3dstuff;

import java.awt.*;
import javax.swing.JPanel;

public class Main extends javax.swing.JFrame {
    
    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    
    static int height = Toolkit.getDefaultToolkit().getScreenSize().height;
    static int width = Toolkit.getDefaultToolkit().getScreenSize().width;
    
    public Main() {
        
        

        

        setTitle("3D stuff");
        setSize(width, height);
        
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());
        
        setContentPane(mainPanel);
        
        mainPanel.add(new MainPanel(), BorderLayout.CENTER);

        
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setVisible(true);
        
    }
    
    public static void main(String[] args) {
        new Main();
    }

}
