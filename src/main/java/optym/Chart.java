package optym;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.JFrame;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.geom.Ellipse2D;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

public class Chart {

    private final JFrame frame;
    private final Object lock = new Object();

    public Chart(String methodName, List<double[]> functionPoints, List<double[]> iterPoints, double[] finalPoint) {

        XYSeriesCollection dataset = getXySeriesCollection(functionPoints, iterPoints, finalPoint);

        JFreeChart chart =  ChartFactory.createXYLineChart(methodName, "x", "f(x)",
                dataset, PlotOrientation.VERTICAL, true, true, false);

        XYPlot plot = chart.getXYPlot();
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();

        renderer.setSeriesShapesVisible(0, false);
        renderer.setSeriesPaint(0, Color.BLUE);
        renderer.setSeriesStroke(0, new BasicStroke(1.0f));

        renderer.setSeriesPaint(1, Color.RED);
        renderer.setSeriesShapesVisible(1, true);
        renderer.setSeriesLinesVisible(1, false);

        renderer.setSeriesPaint(2, Color.MAGENTA);
        renderer.setSeriesShape(2, new Ellipse2D.Double(-5, -5, 10, 10));
        renderer.setSeriesShapesVisible(2, true);
        renderer.setSeriesLinesVisible(2, false);

        plot.setRenderer(renderer);

        ChartPanel chartPanel = new ChartPanel(chart);
        this.frame = new JFrame(methodName);
        this.frame.setContentPane(chartPanel);
        this.frame.setSize(800, 600);
        this.frame.setLocationRelativeTo(null);
        this.frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        this.frame.addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosed(WindowEvent e) {
                synchronized (lock) {
                    lock.notify();
                }
            }
        });
    }

    private static XYSeriesCollection getXySeriesCollection(List<double[]> functionPoints, List<double[]> iterPoints, double[] finalPoint) {
        XYSeries functionSeries = new XYSeries("Function");
        for (double[] pt : functionPoints) {
            functionSeries.add(pt[0], pt[1]);
        }

        XYSeries iterSeries = new XYSeries("Iterations");
        for (double[] pt : iterPoints) {
            iterSeries.add(pt[0], pt[1]);
        }

        XYSeries finalPointSeries = new XYSeries("Ekstremum");
        finalPointSeries.add(finalPoint[0], finalPoint[1]);

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(functionSeries);
        dataset.addSeries(iterSeries);
        dataset.addSeries(finalPointSeries);
        return dataset;
    }

    public void show() {
        this.frame.setVisible(true);
        synchronized (this.lock) {
            try {
                this.lock.wait();
            } catch (InterruptedException ex) {
                ex.printStackTrace();
            }
        }
    }
}
