
    scalar factor = 0.0;

    if (t < 11 && t > 10.0)
    {
        factor = 0.5*(1.0 - cos(2*pi*(t-10.0)));
    }
    if (t < 21 && t > 20.0)
    {
        factor = 0.5*(1.0 - cos(2*pi*(t-20.0)));
    }
    if (t < 31 && t > 30.0)
    {
        factor = 0.5*(1.0 - cos(2*pi*(t-30.0)));
    }
    if (t < 41 && t > 40.0)
    {
        factor = 0.5*(1.0 - cos(2*pi*(t-40.0)));
    }

    if (Pstream::master())
    {
        // rotate moments into body coordinate system
        quaternion rot(theta_.x(), theta_.y(), theta_.z());
        M = rot.invTransform(M);
        
        // apply gravitational force
        F += mass_*g_.value();
        
        // apply spring and damper forces
        F.x() -= (k_.x()*x0_.x()+c_.x()*v0_.x());
        F.y() -= (k_.y()*x0_.y()+c_.y()*v0_.y());
        F.z() -= (k_.z()*x0_.z()+c_.z()*v0_.z());
        M.x() -= (K_.x()*theta0_.x()+C_.x()*omega0_.x());
        M.y() -= (K_.y()*theta0_.y()+C_.y()*omega0_.y());
        M.z() -= (K_.z()*theta0_.z()+C_.z()*omega0_.z());
        
        // apply estimates of added mass force
        F.x() += m_a_.x()*accel_est.x();
        F.y() += m_a_.y()*accel_est.y();
        F.z() += m_a_.z()*accel_est.z();    
        M.x() += I_a_.x()*alpha_est.x();
        M.y() += I_a_.y()*alpha_est.y();
        M.z() += I_a_.z()*alpha_est.z();
        
        // calculate new accelerations
        vector acceleration(vector::zero);
        acceleration.x() = F.x() / ( mass_ + m_a_.x() );
        acceleration.y() = F.y() / ( mass_ + m_a_.y() );
        acceleration.z() = F.z() / ( mass_ + m_a_.z() );
        acceleration.x() *= translationFree_.x();
        acceleration.y() *= translationFree_.y();
        acceleration.z() *= translationFree_.z();
        vector alpha = inv( I_ + symmTensor(I_a_.x(), 0, 0, I_a_.y(), 0, I_a_.z() ) ) & M;
        alpha.x() *= rotationFree_.x();    
        alpha.y() *= rotationFree_.y();    
        alpha.z() *= rotationFree_.z();
        
        vector vTilde, omegaTilde, xTilde, thetaTilde;

        if (order_ <= 1)
        {
            // velocities updated explicitly
            vTilde = v0_ + acceleration*dt_*factor;
            omegaTilde = omega0_ + alpha*dt_*factor;
            
            // positions updated implicity
            xTilde = x0_ + vTilde*dt_*factor;
            thetaTilde = theta0_ + omegaTilde*dt_*factor;
        }
        else
        {
            Info << "WARNING  NOT UPDATED FOR SECOND ORDER!!!!!!" << endl;
            scalar a0 = (dt0_+2.0*dt_)/dt_/(dt0_+dt_);
            scalar a1 = -(dt0_+dt_)/dt0_/dt_;
            scalar a2 = dt_/dt0_/(dt0_+dt_);
            
            // velocities updated explicitly
            vTilde = (-a1*v0_ - a2*v00_ + acceleration)/a0;
            omegaTilde = (-a1*omega0_ - a2*omega00_ + alpha)/a0;
            
            // positions updated implicity
            xTilde = (-a1*x0_ - a2*x00_ + vTilde)/a0;
            thetaTilde = (-a1*theta0_ - a2*theta00_ + omegaTilde)/a0;
        }
        
        // explicitly relax updates
        v_ = beta_*vTilde + (1 - beta_)*v_;
        omega_ = beta_*omegaTilde + (1 - beta_)*omega_;
        x_ = beta_*xTilde + (1 - beta_)*x_;
        theta_ = beta_*thetaTilde + (1 - beta_)*theta_;
        
        // create the recording files if not already done
        makeFiles();
    
        // write out position, velocity, and acceleration
        scalar t = time().value();
        positionFilePtr_() << t << tab << x_ << tab << theta_ << endl;
        velocityFilePtr_() << t << tab << v_ << tab << omega_ << endl;
        accelFilePtr_() << t << tab << (v_-v0_)/dt_ << tab << (omega_-omega0_)/dt_ << endl;
        
        // scatter output to different processors
        for (int ii = 0; ii < Pstream::nProcs(); ii++)
        {
            if (ii != Pstream::myProcNo())
            {
                UOPstream toNeighbor(ii, buf);
                vectorField out(2);
                out[0] = x_;
                out[1] = theta_;
                toNeighbor << out;
            }
        }
        buf.finishedSends();
    }
    else
    {
        buf.finishedSends();
        UIPstream fromMaster(Pstream::masterNo(), buf);
        vectorField in(2, vector::zero);
        
        fromMaster >> in;
        x_ = in[0];
        theta_ = in[1];
    }
