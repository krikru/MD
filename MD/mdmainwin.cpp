////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////

// Standard includes
#include <iostream>

// Own includes
#include "definitions.h"

// Qt includes
#include <QMessageBox>
#include <QCloseEvent>

// Widgets
#include "mdmainwin.h"
#include "ui_mdmainwin.h"

////////////////////////////////////////////////////////////////
// CONSTRUCTOR & DESTRUCTOR
////////////////////////////////////////////////////////////////

mdmainwin::mdmainwin(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::mdmainwin)
{
    // Set up user interface
    ui->setupUi(this);

    // Set up simulation output text browser
    int font_size = 7;
    int tab_width = 8;
    // Create font for the text browser
    QFont log_font("Courier New", font_size + 2, QFont::Normal, false); /* TODO: Why does Qt remove 2 from the font size?? */
    // Create text browser
    QTextEdit *tb = ui->simulation_output_tb;
    tb->setAutoFormatting(QTextEdit::AutoNone); /* Don't format inserted text */
    tb->setFont(log_font);
    tb->setLineWrapMode(QTextEdit::WidgetWidth); /* Wrap text at widget edge */
    tb->setOverwriteMode(false); /* Don't overwrite other text when inserting new one */
    tb->setReadOnly(true); /* Don't accept user inputed text */
    tb->setTabChangesFocus(true); /* Change focus when tab is pressed */
    tb->setTabStopWidth(tab_width * font_size); /* Tab width in pixels */
    tb->setUndoRedoEnabled(false); /* Don't allow undo */
    tb->setWordWrapMode(QTextOption::WrapAnywhere); /* Use all characters places on each line */
}

mdmainwin::~mdmainwin()
{
    delete ui;
}

////////////////////////////////////////////////////////////////
// PRIVATE SLOTS
////////////////////////////////////////////////////////////////

void mdmainwin::on_start_simulation_pb_clicked()
{
#if 1
    // Randomize the simulation with a random seed based on the current time
    uint random_seed = (unsigned int)time(NULL);
#else
    uint random_seed = 0;
#endif
    srand(random_seed);

    // Init element specific constants
#if 0
    //Let's use the Xenon (Xe) atom in an fcc lattice (Melting point 161.4 K)
    //Cohesive energy: 0.16 eV/atom
    //Specific heat: 0.097 J/g

    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(3.98) * P_ANGSTROM;
    ftype epsilon_in = ftype(320e-16) * P_ERG; //1 erg = 10^-7 J
    ftype mass_in = ftype(131.293) * P_U;
    ftype latticeconstant_in = ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 6.200 Å)
    cout<<"Xenon"<<endl;
    // Simulation constants
    ftype dt_in = ftype(0.10) * P_FS; // [s]
    ftype temperature_in = ftype(580.0); // [K]//MSD linear at approx. 800K, why??
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif 1
    //Let's use the Silver (Ag) atom in an fcc lattice (Melting point 1235.08 K) as it is stable at even 500 K
    //Cohesive energy: 2.95 eV/atom

    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(2.65) * P_ANGSTROM;
    ftype epsilon_in = ftype(0.34) * P_EV; //1 erg = 10^-7 J
    ftype mass_in = ftype(107.8682) * P_U;
    ftype latticeconstant_in = ftype(4.090) * P_ANGSTROM;  //  ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 4.090 Å)
    cout<<"Silver"<<endl;
    // Simulation constants
    ftype dt_in = ftype(0.1) * P_FS; // [s]
    ftype temperature_in = ftype(580.0); // [K] MSD linear at approx. 12500 K, why??
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif 0
    //Copper (Melting point 1356.6 K)
    //Cohesive energy: 3.49 eV/atom

    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(2.338) * P_ANGSTROM;
    ftype epsilon_in = ftype(0.4096) * P_EV; //1 erg = 10^-7 J
    ftype mass_in = ftype(63.546) * P_U;
    ftype latticeconstant_in = ftype(3.6100000000) * P_ANGSTROM;//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 3.610 Å)
    cout<<"Copper"<<endl;
    // Simulation constants
    ftype dt_in = ftype(1.0) * P_FS; // [s]
    ftype temperature_in = ftype(800.0); // [K]
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif 1
    //Argon (Melting point 83.8 K)
    //Cohesive energy: 0.080 eV/atom
    //Specific heat: 0.312 J/g

    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(3.40) * P_ANGSTROM;
    ftype epsilon_in = ftype(167e-16) * P_ERG; //1 erg = 10^-7 J
    ftype mass_in = ftype(39.948) * P_U;
    ftype latticeconstant_in = ftype(5.260 * P_ANGSTROM);//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 5.260 Å)
    cout<<"Argon"<<endl;
    // Simulation constants
    ftype dt_in = ftype(0.10) * P_FS; // [s]
    ftype temperature_in = ftype(580.0); // [K]
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#endif

    // Init simulation specific constants
    uint  nrparticles_in = 1000; // The number of particles
    uint  nrinst_in = 1;       // Number of timesteps between measurements of properties
    uint  nrtimesteps_in = 100000; // Desired (or minimum) total number of timesteps
    ftype inner_cutoff_in = ftype(2.5) * sigma_in; //TODO: Make sure this is 2.0 times sigma
    ftype outer_cutoff_in = ftype(1.1) * inner_cutoff_in; //Fewer neighbors -> faster, but too thin skin is not good either. TODO: Change skin thickness to a good one
    uint  impulseresponse_width_in = 10000;                   //decides the number of measured values that are used by the filter function to calculate each output value
    ftype impulseresponse_exponent_in = ftype(0.0);       //the exponent in the impulse response function used to filter the measured values

    // Control
    ftype nrthermostat_time_in = 3;
    ftype thermostat_time_in = nrthermostat_time_in * dt_in;
    ftype deltaEp_in = ftype(0.01);

    // Init flags
    bool thermostat_on_in = !true;
    bool diff_c_on_in = true;
    bool Cv_on_in = true;
    bool pressure_on_in = true;
    bool msd_on_in = true;
    bool Ep_on_in = true;
    bool Ek_on_in = true;

    // Init system and run simulation
    callback<void (*)(void*        )> event_callback_in (static_process_events       , this);
    callback<void (*)(void*, string)> output_callback_in(static_write_to_text_browser, this);
    simulation.set_event_callback (event_callback_in );
    simulation.set_output_callback(output_callback_in);
    simulation.init(nrparticles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_cutoff_in, mass_in, dt_in, nrinst_in, temperature_in, nrtimesteps_in, latticeconstant_in, lattice_type_in, desiredtemp_in, thermostat_time_in, deltaEp_in, impulseresponse_width_in, impulseresponse_exponent_in, thermostat_on_in, diff_c_on_in, Cv_on_in, pressure_on_in, msd_on_in, Ep_on_in, Ek_on_in);
    simulation.run_simulation();
    ui->statusbar->showMessage("Simulation finished.");

    std::cout << "Random seed " << random_seed << std::endl;
}

void mdmainwin::closeEvent(QCloseEvent *event)
{
    if (simulation.is_operating()) {
        // Ask the user whether to abort the operation or not
        QMessageBox msg_box;
        msg_box.setText("An operation is currently being executed.");
        msg_box.setInformativeText("Do you want to abort the operation?");
        msg_box.setStandardButtons(QMessageBox::Abort | QMessageBox::Cancel);
        msg_box.setDefaultButton(QMessageBox::Cancel);
        int result = msg_box.exec();
        if (result != QMessageBox::Abort) {
            event->ignore();
            return;
        }
    }
    simulation.abort_activities();
}

////////////////////////////////////////////////////////////////
// PRIVATE NON-STATIC MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////

void mdmainwin::write_to_text_browser(string output)
{
    QString qstr = QString::fromStdString(output.c_str());

    // Get the pointer to the text browser
    QTextBrowser *tb = ui->simulation_output_tb;
    // Move cursor to insert the text at the end of the text browser
    tb->moveCursor(QTextCursor::End, QTextCursor::MoveAnchor);
    tb->insertPlainText(qstr);
}

////////////////////////////////////////////////////////////////
// PRIVATE STATIC MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////

void mdmainwin::static_write_to_text_browser(void* void_ptr_mainwin, string output)
{
    // Cast to the right pointer type
    mdmainwin *mainwin_ptr = (mdmainwin*)void_ptr_mainwin;

    // Call the non-static version of this function
    mainwin_ptr->write_to_text_browser(output);
}

void mdmainwin::static_process_events(void* void_ptr_mainwin)
{
    // Cast to the right pointer type
    mdmainwin *mainwin_ptr = (mdmainwin*)void_ptr_mainwin;

    // Get the pointer to the system
    mdsystem *sys_ptr = &(mainwin_ptr->simulation);

    if (sys_ptr->is_operating()) {
        uint pre_cent_finished = 100 * sys_ptr->get_loop_num() / sys_ptr->get_max_loops_num();
        mainwin_ptr->ui->statusbar->showMessage("Running simulation... " + QString::number(pre_cent_finished) + " %");
    }
    else {
        mainwin_ptr->ui->statusbar->showMessage("Idle");
    }
    QApplication::processEvents(QEventLoop::AllEvents);
}
