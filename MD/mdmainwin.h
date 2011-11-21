#ifndef MDMAINWIN_H
#define MDMAINWIN_H

#include <QMainWindow>

namespace Ui {
    class mdmainwin;
}

class mdmainwin : public QMainWindow
{
    Q_OBJECT

public:
    explicit mdmainwin(QWidget *parent = 0);
    ~mdmainwin();

private slots:

    void on_start_simulation_pb_clicked();

private:
    // Static member functions
    static void process_events();

private:
    Ui::mdmainwin *ui;
};

#endif // MDMAINWIN_H
