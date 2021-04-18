import React, { Component } from "react";
import Enzyme from "enzyme";
import { connect } from "react-redux";
import Adapter from "enzyme-adapter-react-16";
import { getUsers } from "../../../store/actions/users";

Enzyme.configure({ adapter: new Adapter() });

class Dashboard extends Component {
  render() {
    return (
      <div>
        <section className="content-header">
          <h1>
            Dashboard
            <small>Version 1.0</small>
          </h1>
          <ol className="breadcrumb">
            <li>
              <a href="#">
                <i className="fa fa-dashboard" /> Home
              </a>
            </li>
            <li className="active">Dashboard</li>
          </ol>
        </section>

        <section className="content">
          {/* <!-- Info boxes --> */}
          <div className="row">
            <div className="col-md-3 col-sm-6 col-xs-12">
              <div className="info-box">
                <span className="info-box-icon bg-aqua">
                <i class="fa fa-ellipsis-v" aria-hidden="true"></i>
                </span>

                <div className="info-box-content">
                  <span className="info-box-text">Samples</span>
                  <span className="info-box-number">
                    13
                  </span>
                </div>
                {/* <!-- /.info-box-content --> */}
              </div>
              {/* <!-- /.info-box --> */}
            </div>
            <div className="col-md-3 col-sm-6 col-xs-12">
              <div className="info-box">
                <span className="info-box-icon bg-red">
                  <i className="fa fa-eyedropper" />
                </span>

                <div className="info-box-content">
                  <span className="info-box-text">Tests</span>
                  <span className="info-box-number">1000</span>
                </div>
                {/* <!-- /.info-box-content --> */}
              </div>
              {/* <!-- /.info-box --> */}
            </div>
            {/* <!-- /.col --> */}
            {/* <!-- fix for small devices only --> */}
            <div className="clearfix visible-sm-block" />
            <div className="col-md-3 col-sm-6 col-xs-12">
              <div className="info-box">
                <span className="info-box-icon bg-green">
                  <i className="ion ion-ios-cart-outline" />
                </span>

                <div className="info-box-content">
                  <span className="info-box-text">Sales</span>
                  <span className="info-box-number">760</span>
                </div>
                {/* <!-- /.info-box-content --> */}
              </div>
              {/* <!-- /.info-box --> */}
            </div>
            <div className="col-md-3 col-sm-6 col-xs-12">
              <div className="info-box">
                <span className="info-box-icon bg-yellow">
                  <i className="ion ion-ios-people-outline" />
                </span>

                <div className="info-box-content">
                  <span className="info-box-text">Users</span>
                  <span className="info-box-number">2,000</span>
                </div>
                {/* <!-- /.info-box-content --> */}
              </div>
              {/* <!-- /.info-box --> */}
            </div>
          </div>
        </section>
      </div>
    );
  }
}

export default Dashboard;
