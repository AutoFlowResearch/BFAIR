import React, { Component } from 'react';
import { connect } from 'react-redux';
import { loginUser } from '../../store/actions/auth';

class Login1 extends Component {
  constructor() {
    super();
    this.state = {
      data: {
        usrname: '',
        password: '',
        role: '',
      },
    };
  }

  handleChange = ({ target: input }) => {
    const data = { ...this.state.data };
    data[input.name] = input.value;
    this.setState({ data });
  };

  handleSubmit = async (e) => {
    e.preventDefault();
    const { dispatch } = this.props;
    await dispatch(loginUser('xyz', 'xyz'));
    const { state } = this.props.location;
    window.location = state ? state.from.pathname : '/';
    // this.props.history.push("/");
  };

  render() {
    return (
      <div
        style={{
          position: 'absolute',
          height: '100%',
          width: '100%',
          backgroundPosition: 'center',
          backgroundRepeat: 'no-repeat',
          backgroundSize: 'cover',
          backgroundColor: '#E8E8E8',
        }}
      >
        <div className="login-box">
          <div className="login-logo">
            <a href="../../index2.html">{/* <b>Neo</b>Soft */}</a>
          </div>
          {/* <!-- /.login-logo --> */}
          <div className="login-box-body">
            <p className="login-box-msg">Sign in to start your session</p>

            <form action="../../index2.html" method="post">
              <div className="form-group has-feedback">
                <input
                  type="email"
                  className="form-control"
                  placeholder="Email"
                />
                <span className="glyphicon glyphicon-envelope form-control-feedback" />
              </div>
              <div className="form-group has-feedback">
                <input
                  type="password"
                  className="form-control"
                  placeholder="Password"
                />
                <span className="glyphicon glyphicon-lock form-control-feedback" />
              </div>
              <div className="row">
                <div className="col-xs-8">
                  <div className="checkbox icheck">
                    <label>
                      <input type="checkbox" /> Remember Me
                    </label>
                  </div>
                </div>
                {/* <!-- /.col --> */}
                <div className="col-xs-4">
                  <button
                    type="submit"
                    className="btn btn-primary btn-block btn-flat"
                    onClick={this.handleSubmit}
                  >
                    Sign In
                  </button>
                </div>
                {/* <!-- /.col --> */}
              </div>
            </form>

            <div className="social-auth-links text-center">
              <p>- OR -</p>
              <a
                href="#"
                className="btn btn-block btn-social btn-facebook btn-flat"
              >
                <i className="fa fa-facebook" /> Sign in using Facebook
              </a>
              <a
                href="#"
                className="btn btn-block btn-social btn-google btn-flat"
              >
                <i className="fa fa-google-plus" /> Sign in using Google+
              </a>
            </div>
            {/* <!-- /.social-auth-links --> */}

            <a href="#">I forgot my password</a>
            <br />
            <a href="register.html" className="text-center">
              Register a new membership
            </a>
          </div>
          {/* <!-- /.login-box-body --> */}
        </div>
      </div>
    );
  }
}

export default connect()(Login1);
