import React, { Component } from 'react';
import { connect } from 'react-redux';
import { loginUser } from '../../store/actions/auth';

import Enzyme from 'enzyme';
import Adapter from 'enzyme-adapter-react-16';
import Logo from '../../assets/images/mountain-1.jpg';

Enzyme.configure({ adapter: new Adapter() });

class LoginForm extends Component {
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

  componentDidMount() {
    localStorage.removeItem('user');
  }

  handleSubmit = async (e) => {
    e.preventDefault();
    const { dispatch } = this.props;
    await dispatch(loginUser('xyz', 'xyz'));
    const { state } = this.props.location;
    window.location = state ? state.from.pathname : '/';
    // this.props.history.push("/mainpage");
  };

  render() {
    const { isAuthenticating } = this.props;
    return (
      <div
        style={{
          position: 'absolute',
          height: '100%',
          width: '100%',
          backgroundPosition: 'center',
          backgroundRepeat: 'no-repeat',
          backgroundSize: 'cover',
          backgroundImage: `url(${Logo})`,
        }}
      >
        <div className="col-sm-4">
          <h3>Login</h3>
          <form onSubmit={this.handleSubmit}>
            <div className="form-group">
              <label
                style={{
                  fontSize: 16,
                  fontWeight: 'bold',
                  color: '#fff',
                  opacity: 0.7,
                }}
                htmlFor="username"
              >
                Username
              </label>
              <input
                type="email"
                onChange={this.handleChange}
                className="form-control"
                id="username"
                name="username"
                aria-describedby="emailHelp"
                placeholder="Enter username"
              />
              <small
                style={{
                  fontWeight: 'bold',
                  color: '#fff',
                  opacity: 0.7,
                }}
                id="emailHelp"
                className="form-text text-muted"
              >
                We'll never share your email with anyone else.
              </small>
            </div>
            <div className="form-group">
              <label
                style={{
                  fontSize: 16,
                  fontWeight: 'bold',
                  color: '#fff',
                  opacity: 0.7,
                }}
                htmlFor="password"
              >
                Password
              </label>
              <input
                type="password"
                onChange={this.handleChange}
                className="form-control"
                id="password"
                name="password"
                placeholder="Enter Password"
              />
            </div>
            <div className="form-group">
              <label
                htmlFor="role"
                style={{
                  fontSize: 16,
                  fontWeight: 'bold',
                  color: '#fff',
                  opacity: 0.7,
                }}
              >
                Select Role
              </label>
              <select
                placeholder="select role"
                className="form-control"
                id="role"
                name="role"
                defaultValue="selectrole"
                onChange={this.handleChange}
              >
                <option value="selectrole" disabled>
                  Select Role
                </option>
                <option value="admin">Admin</option>
                <option value="operator">Operator</option>
                <option value="qa">QA</option>
              </select>
            </div>
            <button
              type="submit"
              disabled={isAuthenticating}
              className="btn btn-primary"
            >
              {isAuthenticating ? 'Loading' : 'Submit'}
            </button>
          </form>
        </div>
      </div>
    );
  }
}

const mapStateTOProps = (store) => {
  return { isAuthenticating: store.auth.isAuthenticating };
};

export default connect(mapStateTOProps)(LoginForm);
