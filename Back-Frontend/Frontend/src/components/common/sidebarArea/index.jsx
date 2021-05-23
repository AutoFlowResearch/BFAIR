import React, { Component } from 'react';
import { NavLink } from 'react-router-dom';
import { connect } from 'react-redux';
import { logout } from '../../../store/actions/auth';
import user from '../../../../node_modules/admin-lte/dist/img/user2-160x160.jpg';

class SidebarArea extends Component {
  // handleLogout = () => {
  //   const { dispatch } = this.props;
  //   dispatch(logout());
  //   window.location = "/";
  // };

  render() {
    // const { pathname } = this.props.location;
    // const currentTab = pathname.split("/").pop();
    return (
      <div className="sideBar">
        <nav class="navbar">
          <div className="brand">{/* <h5 className="pl-3">BFAIR</h5>  */}</div>

          <ul className="navbar-nav pt-3">
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Dashboard
              </a>
            </li>
            <li className="nav-item">
              <a class="nav-link" href="#">
                <i class="far fa-file-code"></i>Source code
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Ideas
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Contacts
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Agents
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Articles
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Settings
              </a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="#">
                <i class="far fa-file-code"></i>Subscriptions
              </a>
            </li>
          </ul>
        </nav>
        <div className="sidebarBox">
          <div className="sideboxContent">
            <span>Need Help?</span>
            <p>Please Check our docs</p>
            <button type="button" className="btn leftBtn">
              Documentation
            </button>
          </div>
        </div>
      </div>
    );
  }
}

export default connect()(SidebarArea);

{
  /* <div classNameName="col-sm-3 sidenav">
  <div classNameName="list-group mt-2">
    <NavLink
      to="/mainpage/admin/dashboard"
      classNameName={`list-group-item ${
        currentTab === "dashboard" ? "active" : ""
        }`}
    >
      Dashboard
          </NavLink>

    <NavLink
      to="/mainpage/admin/user-management"
      classNameName={`list-group-item ${
        currentTab === "user-management" ? "active" : ""
        }`}
    >
      User Management
          </NavLink>

    <NavLink
      to="/mainpage/admin/posts"
      classNameName={`list-group-item ${
        currentTab === "posts" ? "active" : ""
        }`}
    >
      Posts
          </NavLink>
    <span classNameName="list-group-item" onClick={this.handleLogout}>
      <button classNameName="btn btn-danger">Logout</button>
    </span>
  </div>
</div> */
}
