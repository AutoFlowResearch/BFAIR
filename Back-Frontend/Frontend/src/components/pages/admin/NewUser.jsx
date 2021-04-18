import React, { Component } from "react";
import { connect } from "react-redux";
import {
  getUserByid,
  updateData,
  clearUserData
} from "../../../store/actions/users";

class UserForm extends Component {
  constructor() {
    super();
    this.handleChange = this.handleChange.bind(this);
  }
  async componentDidMount() {
    const { id } = this.props.match.params;
    if (id === "new-user") return;

    const { dispatch } = this.props;
    await dispatch(getUserByid(id));
  }

  handleChange({ target }) {
    const { dispatch, data } = this.props;
    const { id, value } = target;
    const newdata = { ...data };
    const addressfileds = ["street", "city"];
    if (addressfileds.includes(id)) {
      newdata.address[id] = value;
    } else newdata[id] = value;

    dispatch(updateData(newdata));
  }

  componentWillUnmount() {
    this.props.dispatch(clearUserData());
  }

  handleSubmit = () => {
    console.log("submitted");
  };

  render() {
    const { data, error, isFetching, match } = this.props;
    const { id } = match.params;
    const { name, username, email, address } = data;
    const { street } = address;
    if (error) return <h3>{error.message}</h3>;
    if (isFetching) return <h3>Loading</h3>;
    return (
      <div className="col-sm-6">
        <div className="box box-info">
          <div className="box-header with-border">
            <h3 className="box-title">User Form</h3>
          </div>
          {/* <!-- /.box-header --> */}
          {/* <!-- form start --> */}
          <form className="form-horizontal" onSubmit={this.handleSubmit}>
            <div className="box-body">
              <div className="form-group">
                <label htmlFor="name" className="col-sm-2 control-label">
                  Name
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter name"
                    type="text"
                    className="form-control"
                    id="name"
                    defaultValue={name}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label htmlFor="username" className="col-sm-2 control-label">
                  Username
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter user name"
                    type="text"
                    className="form-control"
                    id="username"
                    defaultValue={username}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label htmlFor="email" className="col-sm-2 control-label">
                  Email
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter email"
                    type="text"
                    className="form-control"
                    id="email"
                    defaultValue={email}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label htmlFor="street" className="col-sm-2 control-label">
                  Street
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter street"
                    type="text"
                    className="form-control"
                    id="street"
                    defaultValue={street}
                    onChange={this.handleChange}
                  />
                </div>
              </div>
            </div>
            <div className="box-footer">
              <button className="btn btn-primary pull-right">
                {id === "new-user" ? "Submit" : "Update"}
              </button>
            </div>
          </form>
        </div>
      </div>
    );
  }
}

const mapStateToProps = state => ({
  isFetching: state.usermanagement.updateuser.isFetching,
  error: state.usermanagement.updateuser.error,
  data: state.usermanagement.updateuser.data
});

export default connect(mapStateToProps)(UserForm);
