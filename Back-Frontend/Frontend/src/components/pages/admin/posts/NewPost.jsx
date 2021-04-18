import React, { Component, Fragment } from "react";
import { connect } from "react-redux";
import { reduxForm, Field } from "redux-form";
import { newField } from "../../../common/newField";

import {
  getPostById,
  updateData,
  clearPostData,
  newPost
} from "../../../../store/actions/posts";

class NewPost extends Component {
  constructor() {
    super();
    this.state = {
      loading: true
    };
  }

  async componentDidMount() {
    const { id } = this.props.match.params;
    if (id === "new-post") return;

    const { getPostById } = this.props;
    await getPostById(id);
  }

  componentWillUnmount() {
    this.props.clearPostData();
  }

  handleChange = ({ target }) => {
    const { updateData, data } = this.props;
    const { id, value } = target;
    const newdata = { ...data };
    newdata[id] = value;
    updateData(newdata);
  };

  handleSubmit = e => {
    e.preventDefault();
    const { title, body } = this.props.data;

    if (title && body) {
      this.props.newPost({ title, body });
      this.props.history.push("/admin/posts");
    }
  };

  render() {
    const { data, error, isFetching, match } = this.props;
    const { id } = match.params;
    const { title, body } = data;
    if (error) return <h3>{error.message}</h3>;
    if (isFetching) return <h3>Loading</h3>;

    return (
      <div className="col-sm-6">
        <div className="box box-info">
          <div className="box-header with-border">
            <h3 className="box-title">Post Form</h3>
          </div>
          <form className="form-horizontal" onSubmit={this.handleSubmit}>
            <div className="box-body">
              <div className="form-group">
                <label htmlFor="title" className="col-sm-2 control-label">
                  Title
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter title"
                    type="text"
                    className="form-control"
                    id="title"
                    defaultValue={title}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label htmlFor="body" className="col-sm-2 control-label">
                  Details
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter details"
                    type="text"
                    className="form-control"
                    id="body"
                    defaultValue={body}
                    onChange={this.handleChange}
                  />
                </div>
              </div>
            </div>
            <div className="box-footer">
              <button
                className="btn btn-primary pull-right"
                disabled={!(title && body)}
              >
                {id === "new-post" ? "Submit" : "Update"}
              </button>
            </div>
          </form>
        </div>
      </div>
    );
  }
}

const mapStateToProps = state => ({
  isFetching: state.postManagement.updatePost.isFetching,
  error: state.postManagement.updatePost.error,
  data: state.postManagement.updatePost.data
});

const mapDispatchToProps = dispatch => ({
  newPost: newpost => dispatch(newPost(newpost)),
  clearPostData: () => dispatch(clearPostData()),
  getPostById: id => dispatch(getPostById(id)),
  updateData: data => dispatch(updateData(data))
});

export default connect(
  mapStateToProps,
  mapDispatchToProps
)(NewPost);
