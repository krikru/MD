#ifndef MDMAINWIN_H
#define MDMAINWIN_H

// Own includes
#include "mdsystem.h"

// Qt includes
#include <QMainWindow>

namespace Ui {
    class mdmainwin;
}

class mdmainwin : public QMainWindow
{
    Q_OBJECT

public:
    // Constructor and destructor
    explicit mdmainwin(QWidget *parent = 0);
    ~mdmainwin();

private slots:

    void on_start_simulation_pb_clicked();
    void closeEvent(QCloseEvent *event);


private:
    // Private functions

    // Static private functions
    static void write_to_text_browser(string output);
    static void process_events();

    // Private variables
    mdsystem simulation;

private:
    Ui::mdmainwin *ui;

};

#endif // MDMAINWIN_H
