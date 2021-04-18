import React, { Component } from "react";
import { connect } from "react-redux";
import { getUsers } from "../../../store/actions/users";

class UserManagement extends Component {
  constructor() {
    super();
    this.state = {};
    this.handleUser = this.handleUser.bind(this);
    this.handleNewUser = this.handleNewUser.bind(this);
  }

  handleUser(id) {
    this.props.history.push(`/admin/user-management/${id}`);
  }

  handleNewUser() {
    this.props.history.push("/admin/user-management/new-user");
  }

  async componentDidMount() {
    const { dispatch } = this.props;
    await dispatch(getUsers());
  }

  render() {
    const { data, error, isFetching } = this.props;
    if (error) return <h3>{error.message}</h3>;
    if (isFetching) return <h3>Loading</h3>;
    return (
      <div>
        <button
          className="btn btn-primary"
          style={{ margin: 9 }}
          onClick={this.handleNewUser}
        >
          New User
        </button>
        <ul className="list-group">
          {data.map(item => (
            <li
              key={item.id}
              className="list-group-item"
              onClick={() => this.handleUser(item.id)}
            >
              {item.name}
            </li>
          ))}
        </ul>
      </div>
    );
  }
}

const mapStateToProps = store => ({
  data: store.usermanagement.userslist.data,
  error: store.usermanagement.userslist.error,
  isFetching: store.usermanagement.userslist.isFetching
});

export default connect(mapStateToProps)(UserManagement);
