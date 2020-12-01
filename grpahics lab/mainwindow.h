#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtGui>
#include <QtCore>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
public slots:
    void Mouse_Pressed();
    void showMousePosition(QPoint& pos);
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_show_axes_clicked();

    void on_Draw_clicked();

    void on_set_point1_clicked();

    void on_set_point2_clicked();

    void on_pushgrid_clicked();

    void on_resetbutton_clicked();

    void draw_dda_algo(QRgb color);

    void on_drawline_clicked();

    void on_simple_line_algorithm_clicked();

    void on_bresenham_clicked();

    void on_drawcirlce_clicked();

    void on_perametric_clicked();

    void on_bresenham_circle_clicked();

    void drawCircle(int xc,int yc, int x1,int y1);

    void ellipsePoint(int centerx,int centery,int x1,int y1);

    void on_ellips_parametric_clicked();

    void on_midpoint_ellips_clicked();

    void floodfill(int fx,int fy,QRgb old_color);

    void on_flood_fill_clicked();

    void initEdgeTable();

    void on_scane_line_clicked();

    void boundaryfill(int fx,int fy,QRgb edge_color,QRgb new_color,QRgb col);

    void on_boundary_clicked();

    void draw_Window();

    void computeCode();

    void cohenSutherlandClip(int x1, int y1,int x2, int y2);

    void on_line_cliping_clicked();

    void on_setcp1_button_clicked();

    void on_setcp2_button_clicked();

    void on_polygonButton_clicked();

    void on_polygon_clicked();

    void on_polygon_draw_clicked();

    void on_translet_button_clicked();

    void on_rotate_button_clicked();

    void on_scaling_clicked();

    void on_reflection_clicked();

    void on_composit_clicked();

    void bezierCurve();

    void on_bezier_curve_clicked();

private:
    Ui::MainWindow *ui;
    QPoint p1,p2;
    QPoint cp1,cp2;
    void point(int,int,QRgb color);
};

#endif // MAINWINDOW_H
