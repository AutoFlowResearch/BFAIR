import React, { Component } from "react";
import { NavLink } from "react-router-dom";
import { connect } from "react-redux";
import { logout } from "../../../store/actions/auth";
import user from "../../../../node_modules/admin-lte/dist/img/user2-160x160.jpg";

class SidebarArea extends Component {
  handleLogout = () => {
    const { dispatch } = this.props;
    dispatch(logout());
    window.location = "/";
  };

  render() {
    const { pathname } = this.props.location;
    const currentTab = pathname.split("/").pop();

    return (
      // < !--Left side column.contains the logo and sidebar-- >
      <aside className="main-sidebar">
        {/* <!-- sidebar: style can be found in sidebar.less --> */}
        <section className="sidebar">
          {/* <!-- Sidebar user panel (optional) --> */}
          {/* <div className="user-panel">
            <div className="pull-left image">
              <img src={user} className="img-circle" alt="User Image" />
            </div>
            <div className="pull-left info">
              <p>User</p>
             
              <a href="#">
                <i className="fa fa-circle text-success" /> Online
              </a>
            </div>
          </div> */}

          {/* <!-- search form (Optional) --> */}
          <form action="#" method="get" className="sidebar-form">
            <div className="input-group">
              <input
                type="text"
                name="q"
                className="form-control"
                placeholder="Search..."
              />
              <span className="input-group-btn">
                <button
                  type="submit" 
                  name="search"
                  id="search-btn"
                  className="btn btn-flat"
                >
                  <i className="fa fa-search" />
                </button>
              </span>
            </div>
          </form>
          {/* <!-- /.search form --> */}

          {/* <!-- Sidebar Menu --> */}
          <ul className="sidebar-menu" data-widget="tree">
            {/* <li className="header">MAIN NAVIGATION</li> */}
            {/* <!-- Optionally, you can add icons to the links --> */}
            <li className={currentTab === "dashboard" ? "active" : ""}>
              <NavLink to="/mainpage/admin/dashboard">
                <i className="fa fa-dashboard" />
                <span>Dashboard</span>
              </NavLink>
            </li>
            <li className={currentTab === "doe-module" ? "active" : ""}>
              <NavLink to="/mainpage/admin/doe-module">
              <i className="fa fa-list" /><span>Doe Module</span>
              </NavLink>
            </li>

            {/* <li className={currentTab === "posts" ? "active" : ""}>
              <NavLink to="/mainpage/admin/posts">
                <i className="fa fa-list" /> <span>Posts</span>
              </NavLink>
            </li>

            <li className={currentTab === "simple-form" ? "active" : ""}>
              <NavLink to="/mainpage/admin/simple-form">
                <i className="fa fa-list" /> <span>Simple form</span>
              </NavLink>
            </li>

            <li className={currentTab === "recipe" ? "active" : ""}>
              <NavLink to="/mainpage/admin/recipe">
                <i className="fa fa-list" /> <span>Recipe</span>
              </NavLink>
            </li> */}
          </ul>
          {/* <!-- /.sidebar-menu --> */}
        </section>
        {/* <!-- /.sidebar --> */}
      </aside>
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
