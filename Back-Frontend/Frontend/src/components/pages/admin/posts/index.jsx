import React, { Component } from 'react';
import { connect } from 'react-redux';
import { getPosts } from '../../../../store/actions/posts';

class Posts extends Component {
  async componentDidMount() {
    await this.props.dispatch(getPosts());
  }

  handleNewPost = () => {
    this.props.history.push('/mainpage/admin/posts/new-post');
  };

  handlePost(id) {
    this.props.history.push(`/mainpage/admin/posts/${id}`);
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
          onClick={this.handleNewPost}
        >
          New Post
        </button>
        <ul className="list-group">
          {data.map((post) => (
            <li
              key={post.id}
              className="list-group-item"
              onClick={() => this.handlePost(post.id)}
            >
              {post.title}
            </li>
          ))}
        </ul>
      </div>
    );
  }
}

const mapStateToProps = ({ postManagement }) => ({
  data: postManagement.postsList.data,
  isFetching: postManagement.postsList.isFetching,
  error: postManagement.postsList.error,
});
export default connect(mapStateToProps)(Posts);
