////////////////////////////////////////////////////////////////
// INCLUDE FILES
////////////////////////////////////////////////////////////////

// Own includes
#include "definitions.h"

// Qt includes
#include <QMessageBox>
#include <QCloseEvent>

// Widgets
#include "mdmainwin.h"
#include "ui_mdmainwin.h"
#include <iostream>

////////////////////////////////////////////////////////////////
// MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////

mdmainwin::mdmainwin(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::mdmainwin)
{
    ui->setupUi(this);
}

mdmainwin::~mdmainwin()
{
    delete ui;
}

void mdmainwin::on_start_simulation_pb_clicked()
{
#if 0
    // Randomize the simulation with a random seed based on the current time
    uint random_seed = (unsigned int)time(NULL);
#else
    uint random_seed = 0;
#endif
    srand(random_seed);

    // Init element specific constants
#if 1
    //Let's use the Xenon (Xe) atom in an fcc lattice (Melting point 161.4 K)
    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(3.98) * P_ANGSTROM;
    ftype epsilon_in = ftype(320e-16) * P_ERG; //1 erg = 10^-7 J
    ftype mass_in = ftype(131.293) * P_U;
    ftype latticeconstant_in = ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);

    // Simulation constants
    ftype dt_in = ftype(1.0) * P_FS; // [s]
    ftype temperature_in = ftype(100.0); // [K]
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif 1
    //Let's use the Silver (Ag) atom in an fcc lattice (Melting point 1235.08 K) as it is stable at even 500 K

    /* Should have the following properties:
     * -------------------------------------
     * Cohesive energy (corrected  ): 55.8 Kcal/mol
     * Cohesive energy (uncurrected): 58.9 Kcal/mol
     * Cohesive energy (observed   ): 68.0 Kcal/mol
     */

    // Element constants
    uint  lattice_type_in = LT_FCC; // (enum_lattice_types)
    ftype sigma_in = ftype(2.65) * P_ANGSTROM;
    ftype epsilon_in = ftype(0.34) * P_EV; //1 erg = 10^-7 J
    ftype mass_in = ftype(107.8682) * P_U;
    ftype latticeconstant_in = ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);

    // Simulation constants
    ftype dt_in = ftype(1.0) * P_FS; // [s]
    ftype temperature_in = ftype(300.0); // [K]
    ftype desiredtemp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#endif

    // Init simulation specific constants
    uint nrparticles_in = 1000; // The number of particles
    uint nrinst_in = 100;       // Number of timesteps between measurements of properties
    uint nrtimesteps_in = 10000; // Desired (or minimum) total number of timesteps
    ftype inner_cutoff_in = ftype(2.0) * sigma_in; //TODO: Make sure this is 2.0 times sigma
    ftype outer_cutoff_in = ftype(1.01) * inner_cutoff_in; //Fewer neighbors -> faster, but too thin skin is not good either. TODO: Change skin thickness to a good one

    // Control
    ftype nrthermostat_time_in = 3;
    ftype thermostat_time_in = nrthermostat_time_in * dt_in;
    ftype deltaEp_in = ftype(0.1);

    // Init flags
    bool thermostat_on_in = !true;
    bool diff_c_on_in = true;
    bool Cv_on_in = true;
    bool pressure_on_in = true;
    bool msd_on_in = true;
    bool Ep_on_in = true;
    bool Ek_on_in = true;

    // Init system and run simulation
    simulation.init(write_to_text_browser, process_events, nrparticles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_cutoff_in, mass_in, dt_in, nrinst_in, temperature_in, nrtimesteps_in, latticeconstant_in, lattice_type_in, desiredtemp_in, thermostat_time_in, deltaEp_in, thermostat_on_in, diff_c_on_in, Cv_on_in, pressure_on_in, msd_on_in, Ep_on_in, Ek_on_in);
    simulation.run_simulation();

    //std::cout << "Random seed " << random_seed << std::endl;
}

void mdmainwin::process_events()
{
    QApplication::processEvents(QEventLoop::AllEvents);
}

void mdmainwin::write_to_text_browser(string output)
{
    QString qstr = QString::fromStdString(output.c_str());
    //ui->simulation_output_tb->append(qstr); //HELP, Cannot use ui in static member, and if I make it non static I get another problem:
                                              //..\mdmainwin.cpp:101: error: no matching function for call to 'mdsystem::init(<unresolved overloaded function type>, void (&)(), uint&,...)'
}                                             // Anyone who have some idea?

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
