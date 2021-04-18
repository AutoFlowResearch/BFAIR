import React, { Component } from "react";
import user from "../../../../node_modules/admin-lte/dist/img/user2-160x160.jpg";
import { logout } from "../../../store/actions/auth";
import { connect } from "react-redux";

class HeaderArea extends Component {
  handleLogout = () => {
    const { dispatch } = this.props;
    dispatch(logout());
    window.location = "/";
  };

  render() {
    return (
      <header className="main-header">
        {/* <!-- Logo --> */}
        <a href="index2.html" className="logo">
          {/* <!-- mini logo for sidebar mini 50x50 pixels --> */}
          <span className="logo-mini">
            <b>A</b>LT
          </span>
          {/* <!-- logo for regular state and mobile devices --> */}
          <span className="logo-lg">
            <b style={{ color: "#fff" }}>Auto</b>
            <font style={{ color: "#fff" }}>Flow</font>
          </span>
        </a>

        {/* <!-- Header Navbar: style can be found in header.less --> */}
        <nav className="navbar navbar-static-top">
          {/* <!-- Sidebar toggle button--> */}
          <a
            href="#"
            className="sidebar-toggle"
            data-toggle="push-menu"
            role="button"
          >
            <span className="sr-only">Toggle navigation</span>
          </a>
          {/* <!-- Navbar Right Menu --> */}
          <div className="navbar-custom-menu">
            <ul className="nav navbar-nav">
              {/* <!-- Messages: style can be found in dropdown.less--> */}
              <li className="dropdown messages-menu">
                <a href="#" className="dropdown-toggle" data-toggle="dropdown">
                  <i className="fa fa-envelope-o" />
                  {/* <span className="label label-success">4</span> */}
                </a>
                <ul className="dropdown-menu">
                  {/* <li className="header">You have 4 messages</li> */}
                </ul>
              </li>
              <li className="dropdown notifications-menu">
                <a href="#" className="dropdown-toggle" data-toggle="dropdown">
                  <i className="fa fa-bell-o" />
                  {/* <span className="label label-warning">10</span> */}
                </a>
              </li>
              <li className="dropdown tasks-menu">
                <a href="#" className="dropdown-toggle" data-toggle="dropdown">
                  <i className="fa fa-flag-o" />
                  {/* <span className="label label-danger">9</span> */}
                </a>
              </li>

              <li className="dropdown user user-menu">
                <a href="#" className="dropdown-toggle" data-toggle="dropdown">
                  <img src={user} className="user-image" alt="User Image" />
                  <span className="hidden-xs">Lawrence</span>
                </a>
                <ul className="dropdown-menu">
                  {/* <!-- User image --> */}
                  <li className="user-header">
                    <img src={user} className="img-circle" alt="User Image" />

                    <p>
                     Lawrence - Developer
                      <small>Member since Nov. 2020</small>
                    </p>
                  </li>
                  {/* <!-- Menu Body --> */}
                  {/* <li className="user-body">
                    <div className="row">
                      <div className="col-xs-4 text-center">
                        <a href="#">Followers</a>
                      </div>
                      <div className="col-xs-4 text-center">
                        <a href="#">Sales</a>
                      </div>
                      <div className="col-xs-4 text-center">
                        <a href="#">Friends</a>
                      </div>
                    </div>
                  </li> */}
                  {/* <!-- Menu Footer--> */}
                  <li className="user-footer">
                    <div className="pull-left">
                      <a href="#" className="btn btn-default btn-flat">
                        Profile
                      </a>
                    </div>
                    <div className="pull-right" onClick={this.handleLogout}>
                      <a href="#" className="btn btn-default btn-flat">
                        Sign out
                      </a>
                    </div>
                  </li>
                </ul>
              </li>
            </ul>
          </div>
        </nav>
      </header>
    );
  }
}

export default connect()(HeaderArea);
