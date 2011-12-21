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
#include <QTimer>

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

    // Start simulation directly when application has finished loading
    QTimer::singleShot(0, this, SLOT(on_start_simulation_pb_clicked()));
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
    if (simulation.is_operating()) {
        // Inform the user that an operation is currently going on
        QMessageBox msg_box;
        msg_box.setText("An operation is currently being executed.");
        msg_box.setInformativeText("Please wait until the current operation has finished.");
        msg_box.setStandardButtons(QMessageBox::Ok);
        msg_box.setDefaultButton(QMessageBox::Ok);
        msg_box.exec();
        return;
    }

    bool very_slowly_changing_total_energy = false;
    // Init element specific constants

    /*
     * Randomization of the simulation (on/off)
     */
#if 1
    // Start unique simulations each time, based on the current time
    uint random_seed = (unsigned int)time(NULL);
#else
    // Start identical simulations each time
    uint random_seed = 0;
#endif
    srand(random_seed);

    /*
     * Select element (xenon, silver, copper or argon)
     */

#define  XENON   1
#define  SILVER  2
#define  COPPER  3
#define  ARGON   4

//#define  ELEMENT  XENON
#define  ELEMENT  SILVER
//#define  ELEMENT  COPPER
//#define  ELEMENT  ARGON

#if  ELEMENT == XENON
    /*********
     * Xenon *
     *********/
    //Let's use the Xenon (Xe) atom in an fcc lattice (Melting point 161.4 K)
    //Cohesive energy: 0.16 eV/atom
    //Specific heat: 0.097 J/(g*K) at 293 K, 0.179 J/(g*K) at 100 K (http://www.springerlink.com/content/p2875753h4661128/fulltext.pdf)

    // Element constants
    ftype sigma_in            = ftype(3.98) * P_SI_ANGSTROM;
    ftype epsilon_in          = ftype(320e-16) * P_SI_ERG;
    ftype mass_in             = ftype(131.293) * P_SI_U;
    uint  lattice_type_in     = LT_FCC; // (enum_lattice_types)
    ftype lattice_constant_in = ftype(6.200 * P_SI_ANGSTROM);//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 6.200 Å)
    cout << "Xenon" << endl;
    // Simulation constants
    ftype temperature_in  = ftype(200.0); // [K]//MSD linear at approx. 800K, why??
    ftype desired_temp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif  ELEMENT == SILVER
    /**********
     * Silver *
     **********/
    //Let's use the Silver (Ag) atom in an fcc lattice (Melting point 1235.08 K) as it is stable at even 500 K
    //Cohesive energy: 2.95 eV/atom
    //Specific heat: 0.233 J/(g*k) at 293 K

    // Element constants
    ftype sigma_in            = ftype(2.65) * P_SI_ANGSTROM;
    ftype epsilon_in          = ftype(0.34) * P_SI_EV;
    ftype mass_in             = ftype(107.8682) * P_SI_U;
    uint  lattice_type_in     = LT_FCC; // (enum_lattice_types)
    ftype lattice_constant_in = ftype(4.090 * P_SI_ANGSTROM);//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 4.090 Å)
    cout << "Silver" << endl;

    // Simulation constants
#if 1 // Keeping a temperature
    ftype temperature_in  = ftype(300.0); // [K] MSD linear at approx. 12500 K, why??
    ftype desired_temp_in = ftype(300.0); // [K]
#elif 1 // Melting
    ftype temperature_in  = ftype(300.0); // [K] MSD linear until approx. 10000 K, why??
    ftype desired_temp_in = ftype(40000); // [K]
    very_slowly_changing_total_energy = true;
#elif 1 // Solidifying
    ftype temperature_in  = ftype(41000.0); // [K]
    ftype desired_temp_in = ftype(-10000); // [K]
    very_slowly_changing_total_energy = true;
#endif
#elif  ELEMENT == COPPER
    /**********
     * Copper *
     **********/
    //Copper (Melting point 1356.6 K)
    //Cohesive energy: 3.49 eV/atom
    //Specific heat: 0.386 J/(g*K) at 293 K

    // Element constants
    ftype sigma_in            = ftype(2.338) * P_SI_ANGSTROM;
    ftype epsilon_in          = ftype(0.4096) * P_SI_EV;
    ftype mass_in             = ftype(63.546) * P_SI_U;
    uint  lattice_type_in     = LT_FCC; // (enum_lattice_types)
    ftype lattice_constant_in = ftype(3.610) * P_SI_ANGSTROM;//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 3.610 Å)
    cout << "Copper" << endl;
    // Simulation constants
    ftype temperature_in  = ftype(580.0); // [K]
    ftype desired_temp_in = temperature_in*ftype(0.9); //TODO: Why times 0.9?
#elif  ELEMENT == ARGON
    /*********
     * Argon *
     *********/
    //Argon (Melting point 83.8 K)
    //Cohesive energy: 0.080 eV/atom
    //Specific heat: 0.312 J/(g*K) at 293 K, approx 0.55 J/(g*K) at 60 K (http://www.springerlink.com/content/k328237200233456/fulltext.pdf)

    // Element constants
    ftype sigma_in            = ftype(3.40) * P_SI_ANGSTROM;
    ftype epsilon_in          = ftype(167e-16) * P_SI_ERG;
    ftype mass_in             = ftype(39.948) * P_SI_U;
    uint  lattice_type_in     = LT_FCC; // (enum_lattice_types)
    ftype lattice_constant_in = ftype(5.260 * P_SI_ANGSTROM);//ftype((pow(2.0, 1.0/6.0)*sigma_in) * M_SQRT2);//(Listed lattice constant 5.260 Å)
    cout << "Argon" << endl;
    // Simulation constants
    ftype temperature_in  = ftype(120.0); // [K]
    ftype desired_temp_in = ftype(100.0); // [K]
#endif

    /*
     * Sampling and filtering of properties
     */
#if  FILTER == KRISTOFERS_FILTER
    uint sample_period_in = 5; // Number of timesteps between each sampling of properties
    ftype default_impulse_response_decay_time_in = ftype(100) * P_SI_FS; //The exponent factor in the impulse response function used to filter the measured values
    /* Select the number of times to apply the filter every time filtering */
    //uint default_num_times_filtering_in = 0; // No filtering.
    uint default_num_times_filtering_in = 1; // Filter once
    //uint default_num_times_filtering_in = 2; // Double filtering
    bool slope_compensate_by_default_in = false;

    uint ensemble_size_in = 0; // Is never used
#elif  FILTER == EMILS_FILTER
    uint sample_period_in = 5;
    uint ensemble_size_in = 100; // Number of values used to calculate averages
    uint default_num_times_filtering_in = 0; // Is never used
    bool slope_compensate_by_default_in = 0; // Is never used

    ftype default_impulse_response_decay_time_in = 0; // Is never used
#endif

    /*
     * Simulation specific constants
     */
    ftype dt_in              = ftype(1.0) * P_SI_FS; // [s]
    uint  num_time_steps_in  = 500; // Desired (or minimum) total number of timesteps
    uint  num_particles_in   = 5000; // The desired (or maximum) number of particles
    ftype thermostat_time_in = ftype(500) * P_SI_FS;
    ftype inner_cutoff_in    = ftype(2.5) * sigma_in; //TODO: Make sure this is 2.0 times sigma
    ftype outer_cutoff_in    = ftype(1.1) * inner_cutoff_in; //Fewer neighbors -> faster, but too thin skin is not good either. TODO: Change skin thickness to a good one
    ftype dEp_tolerance_in   = ftype(1.0);

    /*
     * Simulatin flags
     */
    // Control
    bool thermostat_on_in = true;
    // Measurement
    bool diff_c_on_in   = true;
    bool Cv_on_in       = true;
    bool pressure_on_in = true;
    bool msd_on_in      = true;
    bool Ep_on_in       = true;
    bool Ek_on_in       = true;

    if (very_slowly_changing_total_energy) {
        /*
         * Make sure the final temperature will be approximatelly
         * desired_temp_in no matter how long the thermostat time is.
         */
        thermostat_time_in = 100 * num_time_steps_in * dt_in;
        ftype k = exp(-dt_in*num_time_steps_in/thermostat_time_in);
        desired_temp_in = (desired_temp_in - temperature_in * k)/(1 - k);
    }

    // Init system and run simulation
    callback<void (*)(void*        )> event_callback_in (static_process_events       , this);
    callback<void (*)(void*, string)> output_callback_in(static_write_to_text_browser, this);
    simulation.set_event_callback (event_callback_in );
    simulation.set_output_callback(output_callback_in);
    simulation.init(num_particles_in, sigma_in, epsilon_in, inner_cutoff_in, outer_cutoff_in, mass_in, dt_in, ensemble_size_in, sample_period_in, temperature_in, num_time_steps_in, lattice_constant_in, lattice_type_in, desired_temp_in, thermostat_time_in, dEp_tolerance_in, default_impulse_response_decay_time_in, default_num_times_filtering_in, slope_compensate_by_default_in, thermostat_on_in, diff_c_on_in, Cv_on_in, pressure_on_in, msd_on_in, Ep_on_in, Ek_on_in);
    if (simulation.is_initialized()) {
        simulation.run_simulation();
        ui->statusbar->showMessage("Simulation finished.");
    }
    else {
        ui->statusbar->showMessage("Initialization failed");
    }

    std::cout << "Random seed " << random_seed << std::endl;
    if (ui->close_when_finished_cb->checkState() == Qt::Checked) {
        this->close();
    }
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

void mdmainwin::on_sigma_le_editingFinished()
{
    ui->statusbar->showMessage("Editing sigma finished.");
}

void mdmainwin::on_epsilon_le_editingFinished()
{
    ui->statusbar->showMessage("Editing epsilon finished.");
}

void mdmainwin::on_mass_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_lattice_constant_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_lattice_type_cb_activated(const QString &arg1)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_epsilon_unit_cb_activated(const QString &arg1)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_num_particles_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_init_temperature_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_desired_pressure_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_desire_temperature_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_desired_pressure_unit_cb_activated(const QString &arg1)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_time_step_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_num_time_steps_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_inner_cutoff_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_outer_cutoff_le_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_measurement_interval_sb_editingFinished()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_npe_rb_clicked()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_nve_rb_clicked()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_nvt_rb_clicked()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_npt_rb_clicked()
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_diffoceff_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_pressure_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_cv_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_msd_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_energy_total_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_energy_kinetic_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_energy_potential_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_cohesive_energy_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_store_particle_possitions_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}

void mdmainwin::on_settings_bb_accepted()
{
    ui->statusbar->showMessage("Settings accepted");
}

void mdmainwin::on_settings_bb_rejected()
{
    ui->statusbar->showMessage("Settings rejected");
}

/*
void mdmainwin::on_draw_particles_cb_clicked(bool checked)
{
    ui->statusbar->showMessage("Code needed");
}
*/

void mdmainwin::on_save_element_pb_clicked()
{
    ui->statusbar->showMessage("Saving element... (don't wait in vain)");
}

void mdmainwin::on_load_element_pb_clicked()
{
    ui->statusbar->showMessage("Loading element... (don't wait in vain)");
}
