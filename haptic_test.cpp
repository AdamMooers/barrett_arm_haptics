#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()

#include <iostream>
#include <string>

#include <boost/tuple/tuple.hpp>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>

#include <barrett/standard_main_function.h>


using namespace barrett;
using detail::waitForEnter;

// A System for getting force data from the force and torque sensor
template <size_t DOF>
class ftSystem : public systems::System {
    BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:		Input<double> input;
public:		Output<cf_type> output;
protected:	Output<cf_type>::Value* outputValue;

public:
	explicit ftSystem(ProductManager& pm,
			const std::string& sysName = "ftSystem") :
		    systems::System(sysName), input(this), output(this, &outputValue) {

            fts = NULL;
	        if (pm.foundForceTorqueSensor()) {
		        fts = pm.getForceTorqueSensor();
		        fts->tare();
	        }

        }

	virtual ~ftSystem() { mandatoryCleanUp(); }

protected:
	virtual void operate() {
        fts->update(true);

        //boost::get<0>(ftsIn) = fts->getForce();

        outputValue->setData(&fts->getForce());  // Push data into the output
	}

    ForceTorqueSensor* fts;
    //boost::tuple<cf_type> ftsIn;
};

// Controls the individual joints
template <size_t DOF>
class JpCircle : public systems::SingleIO< typename units::CartesianForce::type, typename units::JointPositions<DOF>::type> {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:
	JpCircle(jp_type startPos, double amplitude, double omega, const std::string& sysName = "JpCircle") :
		systems::SingleIO<cf_type, jp_type>(sysName), jp(startPos), jp_0(jp), amp(amplitude), omega(omega) {
		    i1 = 2;
			i2 = 3;
		}
	virtual ~JpCircle() { this->mandatoryCleanUp(); }

protected:
	jp_type jp;
	jp_type jp_0;
	double amp, omega;
	double theta;
	int i1, i2;

	virtual void operate() {
		theta = this->input.getValue()[2];

		jp[i1] = theta/100.0+jp_0[i1];
		//jp[i2] = amp * (std::cos(theta) - 1.0) + jp_0[i2];

		this->outputValue->setData(&jp);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(JpCircle);
};


template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

    // Configure the force/torque sensor
    ftSystem<DOF> fts(pm);

	wam.gravityCompensate();

	const double TRANSITION_DURATION = 0.5;  // seconds

	const double JP_AMPLITUDE = 0.4;    // radians
	const double CP_AMPLITUDE = 0.1;    // meters
	const double FREQUENCY = 1.0;       // rad/s

	//Rate Limiter
	jp_type rt_jp_cmd;
	systems::RateLimiter<jp_type> jp_rl;

	//Sets the joints to move at 1 m/s
	const double rLimit[] = {1, 1, 1, 1, 1, 1, 1};

	for(size_t i = 0; i < DOF; ++i)
		rt_jp_cmd[i] = rLimit[i];

  // Set start position, depending on robot type and configuration.
	jp_type startPos(0.0);
	if (DOF > 3) {
		// WAM
		startPos[1] = -M_PI_2;
		startPos[3] = M_PI_2 + JP_AMPLITUDE;
	} else {
		std::cout << "Error: No known robot with DOF < 3. Quitting." << std::endl;
		// error
		return -1;
	}

	systems::Ramp time(pm.getExecutionManager(), 1.0);


	printf("Press [Enter] to move the end-point in circles using joint position control.");
	waitForEnter();

	wam.moveTo(startPos);
	//Indicate the current position and the maximum rate limit to the rate limiter
	jp_rl.setCurVal(wam.getJointPositions());
	jp_rl.setLimit(rt_jp_cmd);

	JpCircle<DOF> jpc(startPos, JP_AMPLITUDE, FREQUENCY);

    systems::connect(time.output, fts.input);
    systems::connect(fts.output, jpc.input);

	//Enforces that the individual joints move less than or equal to the above mentioned rate limit
	systems::connect(jpc.output, jp_rl.input);
	wam.trackReferenceSignal(jp_rl.output);
	time.smoothStart(TRANSITION_DURATION);

	printf("Press [Enter] to stop.");
	waitForEnter();
	time.smoothStop(TRANSITION_DURATION);
	wam.idle();

	printf("Press [Enter] to return home.");
	waitForEnter();
	wam.moveHome();
	wam.idle();

	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}
