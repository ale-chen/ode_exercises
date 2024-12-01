package com.ode_exercises;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.event.*;
import java.util.ArrayList;

public class ODEVisualizer {
    // Static wrapper methods for backwards compatibility
    public static void visualize(double[] x, double[] y, double[] z) {
        ArrayList<ArrayList<double[]>> trajectories = new ArrayList<>();
        ArrayList<double[]> singleTrajectory = new ArrayList<>();
        singleTrajectory.add(x);
        singleTrajectory.add(y);
        singleTrajectory.add(z);
        trajectories.add(singleTrajectory);
        MultiParticleVisualizer.visualize(trajectories);
    }
}

class MultiParticleVisualizer extends JFrame {
    public MultiParticleVisualizer(ArrayList<ArrayList<double[]>> trajectories) {
        setTitle("3D ODE Visualization");
        setSize(800, 600);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        add(new Plot3DPanel(trajectories));
    }

    public static void visualize(ArrayList<ArrayList<double[]>> trajectories) {
        SwingUtilities.invokeLater(() -> {
            new MultiParticleVisualizer(trajectories).setVisible(true);
        });
    }

    static class Point3D {
        double x, y, z;
        public Point3D(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    }

    class Plot3DPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {
        private int currentIndex = 0;
        private ArrayList<ArrayList<Point3D>> trajectories = new ArrayList<>();
        private Timer timer;
        private double rotationX = 0;
        private double rotationY = 0;
        private double scale = 200;
        private Point lastMousePos;
        private ArrayList<ArrayList<double[]>> inputTrajectories;
        private boolean isPaused = false;
        private JPanel drawingPanel;
        private double speedMultiplier = 1.0;
        private boolean showAxes = true;
        private boolean showGrid = true;

        public Plot3DPanel(ArrayList<ArrayList<double[]>> inputTrajectories) {
            this.inputTrajectories = inputTrajectories;
            
            for (int i = 0; i < inputTrajectories.size(); i++) {
                trajectories.add(new ArrayList<>());
            }
            
            setLayout(new BorderLayout());

            drawingPanel = new JPanel() {
                @Override
                protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Plot3DPanel.this.paintVisualization(g);
                }
            };
            drawingPanel.addMouseListener(this);
            drawingPanel.addMouseMotionListener(this);
            drawingPanel.addMouseWheelListener(this);
            drawingPanel.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

            JPanel controlPanel = new JPanel();
            JButton resetButton = new JButton("Reset");
            JButton pauseButton = new JButton("Pause");
            JButton toggleAxesButton = new JButton("Toggle Axes");
            JButton toggleGridButton = new JButton("Toggle Grid");
            
            JSlider speedSlider = new JSlider(JSlider.HORIZONTAL, 1, 500, 100);
            speedSlider.addChangeListener(e -> {
                speedMultiplier = speedSlider.getValue() / 100.0;
            });
        
            resetButton.addActionListener(e -> {
                resetVisualization();
                pauseButton.setText("Pause");  // Reset button text too
            });
            
            pauseButton.addActionListener(e -> {
                isPaused = !isPaused;
                pauseButton.setText(isPaused ? "Play" : "Pause");
            });
        
            toggleAxesButton.addActionListener(e -> {
                showAxes = !showAxes;
                drawingPanel.repaint();
            });
        
            toggleGridButton.addActionListener(e -> {
                showGrid = !showGrid;
                drawingPanel.repaint();
            });
            
            JPanel sliderPanel = new JPanel(new BorderLayout());
            sliderPanel.add(new JLabel("Speed: "), BorderLayout.WEST);
            sliderPanel.add(speedSlider, BorderLayout.CENTER);
            
            controlPanel.add(resetButton);
            controlPanel.add(pauseButton);
            controlPanel.add(toggleAxesButton);
            controlPanel.add(toggleGridButton);
            controlPanel.add(sliderPanel);
        
            add(drawingPanel, BorderLayout.CENTER);
            add(controlPanel, BorderLayout.SOUTH);

            timer = new Timer(50, e -> {
                if (!isPaused && currentIndex < getMinTrajectoryLength()) {
                    int steps = Math.max(1, (int)(speedMultiplier));
                    for(int step = 0; step < steps && currentIndex < getMinTrajectoryLength(); step++) {
                        for (int i = 0; i < inputTrajectories.size(); i++) {
                            ArrayList<double[]> currentTraj = inputTrajectories.get(i);
                            trajectories.get(i).add(new Point3D(
                                currentTraj.get(0)[currentIndex],
                                currentTraj.get(1)[currentIndex],
                                currentTraj.get(2)[currentIndex]
                            ));
                        }
                        currentIndex++;
                    }
                    drawingPanel.repaint();
                }
            });
            timer.start();
        }

        private int getMinTrajectoryLength() {
            int minLength = Integer.MAX_VALUE;
            for (ArrayList<double[]> traj : inputTrajectories) {
                minLength = Math.min(minLength, traj.get(0).length);
            }
            return minLength;
        }

        private void resetVisualization() {
            currentIndex = 0;
            for (ArrayList<Point3D> traj : trajectories) {
                traj.clear();
            }
            isPaused = false;
            drawingPanel.repaint();
        }

        public void mousePressed(MouseEvent e) {
            lastMousePos = e.getPoint();
        }

        public void mouseDragged(MouseEvent e) {
            if (lastMousePos != null) {
                int dx = e.getX() - lastMousePos.x;
                int dy = e.getY() - lastMousePos.y;
                
                rotationX += dx * 0.01;
                rotationY += dy * 0.01;
                
                drawingPanel.repaint();
                lastMousePos = e.getPoint();
            }
        }

        @Override
        public void mouseWheelMoved(MouseWheelEvent e) {
            scale *= Math.pow(0.9, e.getWheelRotation());
            drawingPanel.repaint();
        }

        public void mouseReleased(MouseEvent e) { lastMousePos = null; }
        public void mouseClicked(MouseEvent e) {}
        public void mouseEntered(MouseEvent e) {}
        public void mouseExited(MouseEvent e) {}
        public void mouseMoved(MouseEvent e) {}

        private Point3D transform(Point3D p) {
            double x1 = p.x * Math.cos(rotationX) - p.z * Math.sin(rotationX);
            double y1 = p.y;
            double z1 = p.x * Math.sin(rotationX) + p.z * Math.cos(rotationX);

            double x2 = x1;
            double y2 = y1 * Math.cos(rotationY) - z1 * Math.sin(rotationY);
            double z2 = y1 * Math.sin(rotationY) + z1 * Math.cos(rotationY);

            double depth = 4;
            double scaleFactor = depth / (depth + z2);
            
            return new Point3D(
                x2 * scale * scaleFactor,
                y2 * scale * scaleFactor,
                z2
            );
        }

        private void paintVisualization(Graphics g) {
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2.setColor(Color.WHITE);
            g2.fillRect(0, 0, getWidth(), getHeight());

            g2.translate(drawingPanel.getWidth()/2, drawingPanel.getHeight()/2);

            if (showGrid) {
                drawGrid(g2);
            }
            
            if (showAxes) {
                drawAxis(g2, new Point3D(2, 0, 0), "X", Color.RED);
                drawAxis(g2, new Point3D(0, 2, 0), "Y", Color.GREEN);
                drawAxis(g2, new Point3D(0, 0, 2), "Z", Color.BLUE);
            }

            for (int trajIdx = 0; trajIdx < trajectories.size(); trajIdx++) {
                ArrayList<Point3D> trajectory = trajectories.get(trajIdx);
                if (trajectory.size() > 1) {
                    float hue = (float) trajIdx / trajectories.size();
                    
                    for (int i = 1; i < trajectory.size(); i++) {
                        Point3D p1 = transform(trajectory.get(i-1));
                        Point3D p2 = transform(trajectory.get(i));
                        
                        float brightness = 0.5f + ((float)i / trajectory.size()) * 0.5f;
                        g2.setColor(Color.getHSBColor(hue, 1.0f, brightness));
                        
                        g2.draw(new Line2D.Double(p1.x, p1.y, p2.x, p2.y));
                    }
                    
                    Point3D current = transform(trajectory.get(trajectory.size()-1));
                    g2.setColor(Color.getHSBColor(hue, 1.0f, 1.0f));
                    g2.fill(new Ellipse2D.Double(current.x-5, current.y-5, 10, 10));
                }
            }
        }

        private void drawGrid(Graphics2D g2) {
            g2.setColor(new Color(230, 230, 230));
            for (double i = -2; i <= 2; i += 0.5) {
                Point3D p1 = transform(new Point3D(-2, i, 0));
                Point3D p2 = transform(new Point3D(2, i, 0));
                g2.draw(new Line2D.Double(p1.x, p1.y, p2.x, p2.y));
                
                p1 = transform(new Point3D(i, -2, 0));
                p2 = transform(new Point3D(i, 2, 0));
                g2.draw(new Line2D.Double(p1.x, p1.y, p2.x, p2.y));
            }
        }

        private void drawAxis(Graphics2D g2, Point3D end, String label, Color color) {
            Point3D p = transform(end);
            Point3D origin = transform(new Point3D(0, 0, 0));
            g2.setColor(color);
            g2.draw(new Line2D.Double(origin.x, origin.y, p.x, p.y));
            g2.drawString(label, (float)p.x + 10, (float)p.y + 10);
        }
    }
}