#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPixmap>
#include <QImage>
#include <iostream>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintDevice>
#include <QPoint>
#include <math.h>
#include <vector>
#include <complex>

#define maxHt 1800
#define maxWd 1000
#define maxVer 10000

using namespace std;

typedef complex<double> pnt;

double end_point[2][2];

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->x_axis->hide();
    ui->y_axis->hide();
    connect(ui->frame,SIGNAL(Mouse_Pos()),this,SLOT(Mouse_Pressed()));
    connect(ui->frame,SIGNAL(sendMousePosition(QPoint&)),this,SLOT(showMousePosition(QPoint&)));
}

std::vector<std::pair<int,int> > old_vertex;
int clipper_points[4][2];
int mouse_click=0;

MainWindow::~MainWindow()
{
    delete ui;
}
QImage img=QImage(450,450,QImage::Format_RGB888);
int gridsize;
QRgb old_color=qRgb(0,0,0);
std::vector<std::pair<int,int> > vertex_list;

/***************for scane line ************************************************/
typedef struct edgebucket
{
    int ymax;   //max y-coordinate of edge
    float xofymin;  //x-coordinate of lowest edge point updated only in aet
    float slopeinverse;
}EdgeBucket;

typedef struct edgetabletup
{
    // the array will give the scanline number
    // The edge table (ET) with edges entries sorted
    // in increasing y and x of the lower end

    int countEdgeBucket;    //no. of edgebuckets
    EdgeBucket buckets[maxVer];
}EdgeTableTuple;

EdgeTableTuple EdgeTable[maxHt], ActiveEdgeTuple;

void MainWindow::point(int x,int y,QRgb color=qRgb(0,250,250))
{
    int a,b;
    int r=gridsize;
    int a1=x%r;
    int b1=y%r;
    for(a=-a1+1;a<=r-a1-1;a++){
       for(b=-b1+1;b<=r-b1-1;b++){
           img.setPixelColor(x+a,y+b,color);
       }
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

void MainWindow::showMousePosition(QPoint &pos)
{
    ui->mouse_movement->setText(" X : "+QString::number(pos.x()-ui->frame->width()/2)+", Y : "+QString::number(ui->frame->height()/2-pos.y()));
}
void MainWindow::Mouse_Pressed()
{
    ui->mouse_pressed->setText(" X : "+QString::number(ui->frame->x-ui->frame->width()/2)+", Y : "+QString::number(ui->frame->height()/2-ui->frame->y));
    //on_pushgrid_clicked();
    point(ui->frame->x,ui->frame->y);
    ui->x_axis->move(0,ui->frame->y);
    ui->y_axis->move(ui->frame->x,0);
    if(mouse_click==1)
    {
        old_vertex.push_back(std::make_pair(ui->frame->x,ui->frame->y));
    }
}

void MainWindow::on_show_axes_clicked()
{
    if(ui->show_axes->isChecked())
    {
        ui->x_axis->show();
        ui->y_axis->show();
    }
    else{
        ui->x_axis->hide();
        ui->y_axis->hide();
    }
}
void MainWindow::on_set_point1_clicked()
{
    if(ui->draw_line->isChecked()){
        p1.setX(ui->frame->x);
        p1.setY(ui->frame->y);
    }
    end_point[0][0]=ui->frame->x;
    end_point[0][1]=ui->frame->y;
}

void MainWindow::on_set_point2_clicked()
{
    if(ui->draw_line->isChecked()){
        p2.setX(ui->frame->x);
        p2.setY(ui->frame->y);
    }
    end_point[1][0]=ui->frame->x;
    end_point[1][1]=ui->frame->y;
}

void MainWindow::on_Draw_clicked()
{
    int r0=ui->circle_radius->value();
    QPixmap *pix=new QPixmap(450,450);
    QPainter *paint=new QPainter();
    if(ui->draw_circle->isChecked()){
        p1.setX(ui->frame->x);
        p1.setY(ui->frame->y);
        paint->begin(pix);
        paint->setPen(Qt::lightGray);
        paint->drawEllipse(p1,r0,r0);
        ui->frame->setPixmap(*pix);
        paint->end();
    }
    if(ui->draw_line->isChecked()){
        paint->begin(pix);
        paint->setPen(Qt::white);
        paint->drawLine(p1,p2);
        ui->frame->setPixmap(*pix);
        paint->end();
    }
}

void MainWindow::on_resetbutton_clicked()
{
    int a,b;
    for(a=0;a<ui->frame->width()-1;a++){
        for(b=0;b<ui->frame->height();b++){
                img.setPixelColor(a,b,qRgb(0,0,0));
        }
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

void MainWindow::on_pushgrid_clicked()
{
    on_resetbutton_clicked();
    int a,b;
    int k=ui->gridBox->value();
    gridsize=k;
    for(a=0;a<ui->frame->width();a=a+k){
        for(b=0;b<ui->frame->height();b++){
                img.setPixelColor(a,b,qRgb(255,0,0));
        }
    }
    for(a=0;a<ui->frame->height();a++){
        for(b=0;b<ui->frame->width();b=b+k){
                img.setPixelColor(a,b,qRgb(255,0,0));
        }
    }
    int x1=ui->frame->width()/2;
    int y1=ui->frame->height()/2;
    int a1=x1%k;
    int b1=y1%k;
    for(a=-a1+1;a<gridsize-a1;a++){
        for(b=0;b<ui->frame->width();b++){
                img.setPixelColor(b,x1+a,qRgb(255,250,0));
        }
    }
    for(a=-b1+1;a<gridsize-b1;a++){
        for(b=0;b<ui->frame->height();b++){
                img.setPixelColor(y1+a,b,qRgb(255,250,0));
        }
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

/**************************************line drawing algorithms**************************/


//=================================DDA algorithm============================================

void MainWindow::draw_dda_algo(QRgb color=qRgb(0,250,250))
{
    int sourcex=p1.x()/gridsize;
    int sourcey=p1.y()/gridsize;
    int destx=p2.x()/gridsize;
    int desty=p2.y()/gridsize;
    double dx=destx-sourcex;
    double dy=desty-sourcey;
    double steps=fabs(dx)>fabs(dy)?fabs(dx):fabs(dy);

    double nextx=dx/steps;
    double nexty=dy/steps;
    double x=sourcex*gridsize+gridsize/2;
    double y=sourcey*gridsize+gridsize/2;
    for(int i=0;i<=steps;i++)
    {
        point((int)x,(int)y,color);
        x+=nextx*gridsize;
        y+=nexty*gridsize;
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

void MainWindow::on_drawline_clicked()
{
    QPainter painter(&img);
    QPen pen;
    pen.setWidth(1);
    pen.setColor(Qt::yellow);
    draw_dda_algo();
}

//=====================================simple line drawing algorithm==========================================

void MainWindow::on_simple_line_algorithm_clicked()
{
    QPainter painter(&img);
    QPen pen;
    pen.setWidth(1);
    pen.setColor(Qt::yellow);
    int sourcex=p1.x()/gridsize;
    int sourcey=p1.y()/gridsize;
    int destx=p2.x()/gridsize;
    int desty=p2.y()/gridsize;
    float dx=(destx-sourcex);
    float dy=(desty-sourcey);
    if(fabs(dx)>=fabs(dy))
    {
        int y0,x0,xn,x,y;
        float iy;
        y0=(dx>0)?sourcey:desty;
        x0=(dx>0)?sourcex:destx;
        xn=(dx>0)?destx:sourcex;
        float fy=dy/dx;
        iy=y0;
        for(int ix=x0;ix<=xn;ix++)
        {
            x=ix*gridsize+gridsize/2;
            y=iy*gridsize+gridsize/2;
            point(x,y);
            iy+=fy;
        }

    }
    else
    {
        int y0,x0,yn,x,y;
        float ix;
        y0=(dy>0)?sourcey:desty;
        yn=(dy>0)?desty:sourcey;
        x0=(dy>0)?sourcex:destx;
        float fx=dx/dy;
        ix=x0;
        for(int iy=y0;iy<=yn;iy++)
        {
            y=iy*gridsize+gridsize/2;
            x=ix*gridsize+gridsize/2;
            point(x,y,255);
            ix+=fx;
        }
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

//================================bresenham algorithm========================================

void MainWindow::on_bresenham_clicked()
{
    int sourcex=p1.x()/gridsize;
    int sourcey=p1.y()/gridsize;
    int destx=p2.x()/gridsize;
    int desty=p2.y()/gridsize;
    int dx=destx-sourcex;
    int dy=desty-sourcey;

    //int x0=sourcex*gridsize+gridsize/2,y0=sourcey*gridsize+gridsize/2,x1=destx*gridsize+gridsize/2,y1=desty*gridsize+gridsize/2;
    //int dx=p2.x()-p1.x(),dy=p2.y()-p1.y();
    int stepx,stepy;
    if(dy<0){dy=-dy; stepy=-1;} else {stepy=1;}
    if(dx<0){dx=-dx; stepx=-1;} else {stepx=1;}
    if(dx>=dy)
    {
        int x1,p,x,y;
        p=2*dy-dx;
        x1=destx;
        int iy=sourcey;
        int ix=sourcex;
        while(ix!=x1)
        {
            x=ix*gridsize;
            y=iy*gridsize;
            point(x,y);
            if(p<0)
            {
                p=p+2*dy;
            }
            else
            {
                p=p+2*dy-2*dx;
                iy=iy+stepy;
            }
            ix=ix+stepx;
        }
    }
    else
    {
        int y1,p,x,y;
        p=2*dx-dy;
        y1=desty;
        int iy=sourcey;
        int ix=sourcex;
        while(iy!=y1)
        {
            x=ix*gridsize;
            y=iy*gridsize;
            point(x,y);
            if(p<0)
            {
                p+=2*dx;
            }
            else
            {
                p+=2*dx-2*dy;
                ix+=stepx;
            }
            iy+=stepy;
        }
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

/************************************circle drawing algorithms*************************************/

//===========================midpoint algorithm====================================

void MainWindow::on_drawcirlce_clicked()
{
    int radius=ui->circle_radius->value();
        QPainter painter(&img);
        QPen pen;
        pen.setWidth(1);
        pen.setColor(Qt::green);
        painter.setPen(Qt::green);
        int centerx=ui->frame->x,centery=ui->frame->y;
        centerx=(centerx/gridsize)*gridsize+gridsize/2;
        centery=(centery/gridsize)*gridsize+gridsize/2;
        int x1=radius*gridsize,y1=0;
        point(x1+centerx,y1+centery);
        if(radius>0)
        {
            point(-x1+centerx,-y1+centery);
            point(y1+centerx,x1+centery);
            point(-y1+centerx,-x1+centery);
        }

        int p=(1-radius)*gridsize;
        while(x1>y1)
        {
            y1+=gridsize;
            if(p<=0)
                p=p+2*y1+1;
            else
            {
                x1-=gridsize;
                p=p+2*y1-2*x1+1;
            }
            if(x1<y1) break;
            point(x1+centerx,y1+centery);
            point(-x1+centerx,y1+centery);
            point(x1+centerx,-y1+centery);
            point(-x1+centerx,-y1+centery);
            if(x1!=y1)
            {
                point(y1+centerx,x1+centery);
                point(-y1+centerx,x1+centery);
                point(y1+centerx,-x1+centery);
                point(-y1+centerx,-x1+centery);
            }
        }
        ui->frame->setPixmap(QPixmap::fromImage(img));
}

//===================================parametric algorithm=========================================

void MainWindow::on_perametric_clicked()
{
    int x,y,centerx,centery;
    int r=ui->circle_radius->value();
    centerx=ui->frame->x;centery=ui->frame->y;
    r*=gridsize;
    centerx=(centerx/gridsize)*gridsize+gridsize/2;
    centery=(centery/gridsize)*gridsize+gridsize/2;
    for(int i=0;i<360;i++)
    {
        x=centerx+r*cos(float(i*M_PI)/180);
        y=centery+r*sin(float(i*M_PI)/180);
        point(x,y);
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

//=================================bresenham algorithm============================================

void MainWindow::drawCircle(int xc,int yc, int x1,int y1)
{
    QPainter painter(&img);
    QPen pen;
    pen.setWidth(1);
    pen.setColor(Qt::green);
    painter.setPen(Qt::green);

    point(xc+x1, yc+y1);
    point(xc-x1, yc+y1);
    point(xc+x1, yc-y1);
    point(xc-x1, yc-y1);
    point(xc+y1, yc+x1);
    point(xc-y1, yc+x1);
    point(xc+y1, yc-x1);
    point(xc-y1, yc-x1);
}

void MainWindow::on_bresenham_circle_clicked()
{
    int radius=ui->circle_radius->value();
        QPainter painter(&img);
        QPen pen;
        pen.setWidth(1);
        pen.setColor(Qt::green);
        painter.setPen(Qt::green);
        int centerx=ui->frame->x,centery=ui->frame->y;
        centerx=(centerx/gridsize)*gridsize+gridsize/2;
        centery=(centery/gridsize)*gridsize+gridsize/2;

        int x1 = 0, y1 = radius*gridsize;
        int decision_parameter = (3 - 2 * radius)*gridsize;
        while (y1 >= x1)
        {
            drawCircle(centerx,centery, x1, y1);
            x1+=gridsize;

            if (decision_parameter > 0)
            {
                y1-=gridsize;
                decision_parameter = decision_parameter + 4 * (x1 - y1) + 10;
            }
            else
                decision_parameter = decision_parameter + 4 * x1 + 6;
            drawCircle(centerx,centery, x1, y1);
        }

        ui->frame->setPixmap(QPixmap::fromImage(img));
}

/**********************************************ellips drawing algorithms**********************************/


//===================================parametric algorithm================================================

void MainWindow::on_ellips_parametric_clicked()
{
    int x,y,centerx,centery,rx,ry;
    rx=ui->el_rad1->value();
    ry=ui->el_rad2->value();
    centerx=ui->frame->x;centery=ui->frame->y;
    rx*=gridsize;
    ry*=gridsize;
    centerx=(centerx/gridsize)*gridsize+gridsize/2;
    centery=(centery/gridsize)*gridsize+gridsize/2;
    for(int i=0;i<360;i++)
    {
        x=centerx+rx*cos(float(i*M_PI)/180);
        y=centery+ry*sin(float(i*M_PI)/180);
        point(x,y,255);
    }
    ui->frame->setPixmap(QPixmap::fromImage(img));
}

//==================================midpoint algorithm================================

void MainWindow::ellipsePoint(int centerx,int centery,int x1,int y1)
{
    point(centerx+x1*gridsize,centery+y1*gridsize,255);
    point(centerx+x1*gridsize,centery-y1*gridsize,255);
    point(centerx-x1*gridsize,centery+y1*gridsize,255);
    point(centerx-x1*gridsize,centery-y1*gridsize,255);
}
void MainWindow::on_midpoint_ellips_clicked()
{
    int x1, y1, p;

    int rx=ui->el_rad1->value();
    int ry=ui->el_rad2->value();

    int xc=ui->frame->x,yc=ui->frame->y;
    xc=(xc/gridsize)*gridsize+gridsize/2;
    yc=(yc/gridsize)*gridsize+gridsize/2;

    x1=0;
    y1=ry;
    p=(ry*ry)-(rx*rx*ry)+((rx*rx)/4);
    while((2*x1*ry*ry)<(2*y1*rx*rx))
    {
        ellipsePoint(xc,yc,x1,y1);

        if(p<0)
        {
            x1=x1+1;
            p=p+(2*ry*ry*x1)+(ry*ry);
        }
        else
        {
            x1=x1+1;
            y1=y1-1;
            p=p+(2*ry*ry*x1+ry*ry)-(2*rx*rx*y1);
        }
    }
    p=((float)x1+0.5)*((float)x1+0.5)*ry*ry+(y1-1)*(y1-1)*rx*rx-rx*rx*ry*ry;

    while(y1>=0)
    {
        ellipsePoint(xc,yc,x1,y1);

        if(p>0)
        {
            y1=y1-1;
            p=p-(2*rx*rx*y1)+(rx*rx);
        }
        else
        {
            y1=y1-1;
            x1=x1+1;
        p=p+(2*ry*ry*x1)-(2*rx*rx*y1)-(rx*rx);
        }
   }
}

//=========================================filling============================================================

/**********************************************flood fill**************************************************/

void MainWindow::floodfill(int fx,int fy,QRgb old_color)
{
    if(img.pixel(fx,fy)==old_color || img.pixel(fx,fy)==qRgb(255,250,0))
    {
        point(fx,fy);
        floodfill(fx+gridsize,fy,old_color);
        floodfill(fx-gridsize,fy,old_color);
        floodfill(fx,fy+gridsize,old_color);
        floodfill(fx,fy-gridsize,old_color);
        //floodfill(fx+gridsize,fy-gridsize,old_color);
        //floodfill(fx-gridsize,fy+gridsize,old_color);
        //floodfill(fx-gridsize,fy-gridsize,old_color);
        //floodfill(fx+gridsize,fy+gridsize,old_color);
    }
}
void MainWindow::on_flood_fill_clicked()
{
    int fx=ui->frame->x,fy=ui->frame->y;
    point(fx,fy,old_color);
    ui->frame->setPixmap(QPixmap::fromImage(img));
    fx=(fx/gridsize)*gridsize+gridsize/2;
    fy=(fy/gridsize)*gridsize+gridsize/2;
    floodfill(fx,fy,old_color);
}

/**************************************scane line fill***********************************************/

void MainWindow::initEdgeTable()
{
    int i;
    for (i=0; i<maxHt; i++)
    {
        EdgeTable[i].countEdgeBucket = 0;
    }

    ActiveEdgeTuple.countEdgeBucket = 0;
}


void insertionSort(EdgeTableTuple *ett)
{
    int i,j;
    EdgeBucket temp;

    for (i = 1; i < ett->countEdgeBucket; i++)
    {
        temp.ymax = ett->buckets[i].ymax;
        temp.xofymin = ett->buckets[i].xofymin;
        temp.slopeinverse = ett->buckets[i].slopeinverse;
        j = i - 1;

    while ((temp.xofymin < ett->buckets[j].xofymin) && (j >= 0))
    {
        ett->buckets[j + 1].ymax = ett->buckets[j].ymax;
        ett->buckets[j + 1].xofymin = ett->buckets[j].xofymin;
        ett->buckets[j + 1].slopeinverse = ett->buckets[j].slopeinverse;
        j = j - 1;
    }
    ett->buckets[j + 1].ymax = temp.ymax;
    ett->buckets[j + 1].xofymin = temp.xofymin;
    ett->buckets[j + 1].slopeinverse = temp.slopeinverse;
    }
}
void storeEdgeInTuple (EdgeTableTuple *receiver,int ym,int xm,float slopInv)
{
    (receiver->buckets[(receiver)->countEdgeBucket]).ymax = ym;
    (receiver->buckets[(receiver)->countEdgeBucket]).xofymin = (float)xm;
    (receiver->buckets[(receiver)->countEdgeBucket]).slopeinverse = slopInv;

    insertionSort(receiver);

    (receiver->countEdgeBucket)++;


}

void removeEdgeByYmax(EdgeTableTuple *Tup,int yy)
{
    int i,j;
    for (i=0; i< Tup->countEdgeBucket; i++)
    {
        if (Tup->buckets[i].ymax == yy)
        {
            for ( j = i ; j < Tup->countEdgeBucket -1 ; j++ )
                {
                Tup->buckets[j].ymax =Tup->buckets[j+1].ymax;
                Tup->buckets[j].xofymin =Tup->buckets[j+1].xofymin;
                Tup->buckets[j].slopeinverse = Tup->buckets[j+1].slopeinverse;
                }
                Tup->countEdgeBucket--;
            i--;
        }
    }
}

void updatexbyslopeinv(EdgeTableTuple *Tup)
{
    int i;

    for (i=0; i<Tup->countEdgeBucket; i++)
    {
        (Tup->buckets[i]).xofymin =(Tup->buckets[i]).xofymin + (Tup->buckets[i]).slopeinverse;
    }
}

void MainWindow::on_scane_line_clicked()
{
    int i, j, x1, ymax1, x2, ymax2, FillFlag = 0, coordCount;

    for (i=0; i<maxHt; i++)
    {
        for (j=0; j<EdgeTable[i].countEdgeBucket; j++)
        {
            storeEdgeInTuple(&ActiveEdgeTuple,EdgeTable[i].buckets[j].
                     ymax,EdgeTable[i].buckets[j].xofymin,
                    EdgeTable[i].buckets[j].slopeinverse);
        }

        removeEdgeByYmax(&ActiveEdgeTuple, i);

        insertionSort(&ActiveEdgeTuple);

        j = 0;
        FillFlag = 0;
        coordCount = 0;
        x1 = 0;
        x2 = 0;
        ymax1 = 0;
        ymax2 = 0;
        while (j<ActiveEdgeTuple.countEdgeBucket)
        {
            if (coordCount%2==0)
            {
                x1 = (int)(ActiveEdgeTuple.buckets[j].xofymin);
                ymax1 = ActiveEdgeTuple.buckets[j].ymax;
                if (x1==x2)
                {
                    if (((x1==ymax1)&&(x2!=ymax2))||((x1!=ymax1)&&(x2==ymax2)))
                    {
                        x2 = x1;
                        ymax2 = ymax1;
                    }

                    else
                    {
                        coordCount++;
                    }
                }

                else
                {
                        coordCount++;
                }
            }
            else
            {
                x2 = (int)ActiveEdgeTuple.buckets[j].xofymin;
                ymax2 = ActiveEdgeTuple.buckets[j].ymax;

                FillFlag = 0;
                if (x1==x2)
                {
                    if (((x1==ymax1)&&(x2!=ymax2))||((x1!=ymax1)&&(x2==ymax2)))
                    {
                        x1 = x2;
                        ymax1 = ymax2;
                    }
                    else
                    {
                        coordCount++;
                        FillFlag = 1;
                    }
                }
                else
                {
                    coordCount++;
                    FillFlag = 1;
                }

            if(FillFlag)
            {
                    p1.setX(x1);p1.setY(i);
                    p2.setX(x2);p2.setY(i);
                    on_drawline_clicked();
            }

        }

        j++;
    }
    updatexbyslopeinv(&ActiveEdgeTuple);
}

    vertex_list.clear();
    initEdgeTable();
}

/********************************************boundary fill****************************************/

void MainWindow::boundaryfill(int fx,int fy,QRgb edge_color,QRgb new_color,QRgb col)
{
    if(img.pixel(fx,fy)!=edge_color&&img.pixel(fx,fy)!=new_color)
    {
        point(fx,fy,col);
        boundaryfill(fx+gridsize,fy,edge_color,new_color,col);
        boundaryfill(fx-gridsize,fy,edge_color,new_color,col);
        boundaryfill(fx,fy+gridsize,edge_color,new_color,col);
        boundaryfill(fx,fy-gridsize,edge_color,new_color,col);
        ui->frame->setPixmap(QPixmap::fromImage(img));
    }
}

void MainWindow::on_boundary_clicked()
{
    int fx=ui->frame->x,fy=ui->frame->y;
    point(fx,fy,old_color);
    QRgb edge_color=qRgba(0,250,250,0xFF);
    QRgb new_color=qRgba(250,250,250,0xFF);
    fx=(fx/gridsize)*gridsize+gridsize/2;
    fy=(fy/gridsize)*gridsize+gridsize/2;
    boundaryfill(fx,fy,edge_color,new_color,qRgb(250,250,250));
}

//==========================================Line Cliping==============================================================

void MainWindow::on_setcp1_button_clicked()
{
    cp1.setX((ui->frame->x/gridsize)*gridsize+gridsize/2);
    cp1.setY((ui->frame->y/gridsize)*gridsize+gridsize/2);
}

void MainWindow::on_setcp2_button_clicked()
{
    cp2.setX((ui->frame->x/gridsize)*gridsize+gridsize/2);
    cp2.setY((ui->frame->y/gridsize)*gridsize+gridsize/2);

    clipper_points[0][0]=cp1.x();
    clipper_points[0][1]=cp1.y();
    clipper_points[1][0]=cp1.x();
    clipper_points[1][1]=cp2.y();
    clipper_points[2][0]=cp2.x();
    clipper_points[2][1]=cp2.y();
    clipper_points[3][0]=cp2.x();
    clipper_points[3][1]=cp1.y();

    draw_Window();
}


const int INSIDE = 0; // 0000
const int LEFT = 1;   // 0001
const int RIGHT = 2;  // 0010
const int BOTTOM = 4; // 0100
const int TOP = 8;    // 1000

void MainWindow:: draw_Window()
{
    p1.setX(clipper_points[0][0]);
    p1.setY(clipper_points[0][1]);
    p2.setX(clipper_points[1][0]);
    p2.setY(clipper_points[1][1]);
    draw_dda_algo(qRgb(0,250,250));

    p1.setX(clipper_points[1][0]);
    p1.setY(clipper_points[1][1]);
    p2.setX(clipper_points[2][0]);
    p2.setY(clipper_points[2][1]);
    draw_dda_algo(qRgb(0,250,250));

    p1.setX(clipper_points[2][0]);
    p1.setY(clipper_points[2][1]);
    p2.setX(clipper_points[3][0]);
    p2.setY(clipper_points[3][1]);
    draw_dda_algo(qRgb(0,250,250));

    p1.setX(clipper_points[3][0]);
    p1.setY(clipper_points[3][1]);
    p2.setX(clipper_points[0][0]);
    p2.setY(clipper_points[0][1]);
    draw_dda_algo(qRgb(0,250,250));
}
// Function to compute region code for a point(x, y)

double x_max,y_max,x_min,y_min;
void MainWindow::computeCode()
{
    int x1=cp1.x(),x2=cp2.x(),y1=cp1.y(),y2=cp2.y();
    x_max=(x1>x2)?x1:x2;
    y_max=(y1>y2)?y1:y2;
    x_min=(x1<x2)?x1:x2;
    y_min=(y1<y2)?y1:y2;
}
// Implementing Cohen-Sutherland algorithm
// Clipping a line from P1 = (x2, y2) to P2 = (x2, y2)
void MainWindow::cohenSutherlandClip(int ix1, int iy1,int ix2, int iy2)
{
    int code1=INSIDE,code2=INSIDE;
    computeCode();
    if(ix1<x_min)code1|=LEFT;
    else if(ix1>x_max)code1|=RIGHT;

    if(iy1>y_max)code1|=BOTTOM;
    else if(iy1<y_min)code1|=TOP;

    if(ix2<x_min)code2|=LEFT;
    else if(ix2>x_max)code2|=RIGHT;

    if(iy2>y_max)code2|=BOTTOM;
    else if(iy2<y_min)code2|=TOP;

    std::cout<<code1<<" "<<code2;

    double x1,x2,y1,y2;
    x1=ix1,x2=ix2,y1=iy1,y2=iy2;
    bool accept = false;

        while(true)
        {
            if ((code1 == 0) && (code2 == 0))
            {
                // If both endpoints lie within rectangle
                accept = true;
                break;
            }
            else if (code1 & code2)
            {
                // If both endpoints are outside rectangle,
                // in same region
                break;
            }
            else
            {
                // Some segment of line lies within the
                // rectangle
                int code_out;
                double x, y;
                // At least one endpoint is outside the
                // rectangle, pick it.
                if (code1 != 0)
                    code_out = code1;
                else
                    code_out = code2;
                // Find intersection point;
                // using formulas y = y1 + slope * (x - x1),
                // x = x1 + (1 / slope) * (y - y1)
                if (code_out & TOP)
                {
                    // point is above the clip rectangle
                    x = x1 + (x2 - x1) * (y_min - y1) / (y2 - y1);
                    y = y_min;
                }
                else if (code_out & BOTTOM)
                {
                    // point is below the rectangle
                    x = x1 + (x2 - x1) * (y_max - y1) / (y2 - y1);
                    y = y_max;
                }
                else if (code_out & RIGHT)
                {
                    // point is to the right of rectangle
                    y = y1 + (y2 - y1) * (x_max - x1) / (x2 - x1);
                    x = x_max;
                }
                else if (code_out & LEFT)
                {
                    // point is to the left of rectangle
                    y = y1 + (y2 - y1) * (x_min - x1) / (x2 - x1);
                    x = x_min;
                }

                if (code_out == code1)
                {
                    x1 = x;
                    y1 = y;
                    code1=INSIDE;
                    if(x1<x_min)code1|=LEFT;
                    if(x1>x_max)code1|=RIGHT;
                    if(y1>y_max)code1|=BOTTOM;
                    if(y1<y_min)code1|=TOP;
                }
                else
                {
                    x2 = x;
                    y2 = y;
                    code2=INSIDE;
                    if(x2<x_min)code2|=LEFT;
                    if(x2>x_max)code2|=RIGHT;
                    if(y2>y_max)code2|=BOTTOM;
                    if(y2<y_min)code2|=TOP;
                }
            }
        }
        if (accept)
        {
            //If accepted
            //Just reset and draw the boundary and the line
            //Reset the screen and draw the grid

            cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;

            p1.setX(x1);
            p1.setY(y1);

            p2.setX(x2);
            p2.setY(y2);

            draw_dda_algo(qRgb(0,250,250));
            draw_Window();
    }
    else
    {
        //If not accepted
        //Just reset and draw the boundary
        //Reset the screen and draw the grid
        cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        draw_Window();
    }

}
void MainWindow::on_line_cliping_clicked()
{
    draw_dda_algo(qRgb(0,0,0));
    cohenSutherlandClip(p1.x(),p1.y(),p2.x(),p2.y());
}


/************************************************polygon clipping********************************************/
void MainWindow::on_polygonButton_clicked()
{
    int x,n;
    n=old_vertex.size();
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
        cohenSutherlandClip(p1.x(),p1.y(),p2.x(),p2.y());
    }
}

void MainWindow::on_polygon_clicked()
{
    old_vertex.clear();
    mouse_click=1;
}

void MainWindow::on_polygon_draw_clicked()
{
    int x,n;
    n=old_vertex.size();
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,250,250));
    }
    mouse_click=0;
}


/*********************************************translation,scaling,rotation,reflection,composite Transformation*********************/

//===============================================translation===============================================================


void MainWindow::on_translet_button_clicked()
{
    int k=gridsize;
    int tx=ui->transletx->value();
    int ty=ui->translety->value();
    tx*=k;
    ty*=k;
    ty=-ty;
    int x,n;
    n=old_vertex.size();
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
    }
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first+tx);
        p1.setY(old_vertex[i].second+ty);
        p2.setX(old_vertex[x].first+tx);
        p2.setY(old_vertex[x].second+ty);
        draw_dda_algo(qRgb(250,250,250));
    }
}

//====================================================rotattion=======================================================

void MainWindow::on_rotate_button_clicked()
{
    int angl=ui->angle->value();
    int x,n;
    n=old_vertex.size();
    double dang=(double)angl*M_PI/180.0;
    double sinang=sin(dang);
    double cosang=cos(dang);
    int x1,y1;
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
    }
    int f_org_x=old_vertex[0].first - ui->frame->width()/2;
    int f_org_y=ui->frame->width()/2 - old_vertex[0].second;
    for(int i=0;i<n;i++)
    {
        old_vertex[i].first  =old_vertex[i].first - ui->frame->width()/2-f_org_x;
        old_vertex[i].second = ui->frame->width()/2 - old_vertex[i].second-f_org_y;
    }
    for(int i=0;i<n;i++)
    {
        x1=old_vertex[i].first*cosang - old_vertex[i].second*sinang;
        y1=old_vertex[i].first*sinang + old_vertex[i].second*cosang;

        x1=x1+ ui->frame->width()/2+f_org_x;
        y1= ui->frame->width()/2 - y1-f_org_y;

        old_vertex[i].first  = x1;
        old_vertex[i].second = y1;
    }
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(250,250,250));
    }
}

//===============================================Scaling=================================================

void MainWindow::on_scaling_clicked()
{
    double tx=ui->sx->value();
    double ty=ui->sy->value();

    //int k=gridsize;
    int x,n;
    n=old_vertex.size();
    int x1,y1;
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
    }
    for(int i=0;i<n;i++)
    {
        old_vertex[i].first  =old_vertex[i].first - ui->frame->width()/2;
        old_vertex[i].second = ui->frame->width()/2 - old_vertex[i].second;
    }
    for(int i=0;i<n;i++)
    {
        //x=(i+1)%n;
        x1=old_vertex[i].first*tx;
        y1=old_vertex[i].second*ty;

        x1=x1+ ui->frame->width()/2;
        y1= ui->frame->width()/2 - y1;

        old_vertex[i].first  = x1;
        old_vertex[i].second = y1;
    }
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(250,250,250));
    }
}

//===============================================reflection===========================================

pnt reflect(int a,int b,pnt A,pnt B)
{
   pnt P(a,b);
   pnt pr=P-A;
   pnt br=B-A;
   pnt ref=pr/br;
   return conj(ref)*br+A;
}

void MainWindow::on_reflection_clicked()
{
    pnt A(end_point[0][0],end_point[0][1]);
    pnt B(end_point[1][0],end_point[1][1]);
    pnt P;
    int x,n;
    n=old_vertex.size();
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
    }
    for(int i=0;i<n;i++)
    {
        P=reflect(old_vertex[i].first,old_vertex[i].second,A,B);
        old_vertex[i].first=P.real();
        old_vertex[i].second=P.imag();
    }
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(250,250,250));
    }
}

//=======================================================composit================================================



void MainWindow::on_composit_clicked()
{
    double tx=ui->sx->value();
    double ty=ui->sy->value();
    int angl=ui->angle->value();
    int t1x=ui->transletx->value();
    int t1y=ui->translety->value();
    int k=gridsize;
    int x1,y1;
    t1x*=k;
    t1y*=k;
    t1y=-t1y;
    pnt A(end_point[0][0],end_point[0][1]);
    pnt B(end_point[1][0],end_point[1][1]);
    pnt P;
    int x,n;
    double dang=(double)angl*M_PI/180.0;
    double sinang=sin(dang);
    double cosang=cos(dang);
    n=old_vertex.size();
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first);
        p1.setY(old_vertex[i].second);
        p2.setX(old_vertex[x].first);
        p2.setY(old_vertex[x].second);
        draw_dda_algo(qRgb(0,0,0));
    }
    for(int i=0;i<n;i++)
    {
        P=reflect(old_vertex[i].first,old_vertex[i].second,A,B);
        old_vertex[i].first=P.real();
        old_vertex[i].second=P.imag();
    }
    for(int i=0;i<n;i++)
    {
        //x=(i+1)%n;
        old_vertex[i].first  =old_vertex[i].first - ui->frame->width()/2;
        old_vertex[i].second = ui->frame->width()/2 - old_vertex[i].second;
        x1=old_vertex[i].first*tx;
        y1=old_vertex[i].second*ty;

        x1=x1+ ui->frame->width()/2;
        y1= ui->frame->width()/2 - y1;

        old_vertex[i].first  = x1;
        old_vertex[i].second = y1;
    }
    int f_org_x=old_vertex[0].first - ui->frame->width()/2;
    int f_org_y=ui->frame->width()/2 - old_vertex[0].second;
    for(int i=0;i<n;i++)
    {
        old_vertex[i].first  =old_vertex[i].first - ui->frame->width()/2-f_org_x;
        old_vertex[i].second = ui->frame->width()/2 - old_vertex[i].second-f_org_y;

        x1=old_vertex[i].first*cosang - old_vertex[i].second*sinang;
        y1=old_vertex[i].first*sinang + old_vertex[i].second*cosang;

        x1=x1+ ui->frame->width()/2+f_org_x;
        y1= ui->frame->width()/2 - y1-f_org_y;

        old_vertex[i].first  = x1;
        old_vertex[i].second = y1;
    }
    for(int i=0;i<n;i++)
    {
        x=(i+1)%n;
        p1.setX(old_vertex[i].first+t1x);
        p1.setY(old_vertex[i].second+t1y);
        p2.setX(old_vertex[x].first+t1x);
        p2.setY(old_vertex[x].second+t1y);
        draw_dda_algo(qRgb(250,250,250));
    }
}

/******************************************************Bezier curve *****************************************************/


void MainWindow::bezierCurve()
{
    double xu = 0.0 , yu = 0.0 , u = 0.0 ;
    //int i = 0 ;
    for(u = 0.0 ; u <= 1.0 ; u += 0.0001)
    {
        xu = pow(1-u,3)*old_vertex[0].first+3*u*pow(1-u,2)*old_vertex[1].first+3*pow(u,2)*(1-u)*old_vertex[2].first+pow(u,3)*old_vertex[3].first;
        yu = pow(1-u,3)*old_vertex[0].second+3*u*pow(1-u,2)*old_vertex[1].second+3*pow(u,2)*(1-u)*old_vertex[2].second+pow(u,3)*old_vertex[3].second;
        point((int)xu , (int)yu,qRgb(450,0,0));
    }
}

void MainWindow::on_bezier_curve_clicked()
{
   bezierCurve();
}
